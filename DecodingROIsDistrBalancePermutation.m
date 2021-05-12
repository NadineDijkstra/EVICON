function DecodingROIsDistrBalancePermutation(cfg)
% function DecodingROIsDistrBalancePermutation(cfg)

% outputDir
outputDir = fullfile(cfg.root,'Results','GroupResults',[cfg.outputDir '/Permutation']);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% get the ROI indices
load(fullfile(cfg.root,cfg.ROIs),'idx','names');

% loop over subjects
nsubjects  = length(cfg.subjects);
nROIs      = length(names);
permAcc    = zeros(nsubjects,nROIs,cfg.nPerm,cfg.nPerm);
corrY      = cell(nsubjects,cfg.nPerm);
ridx       = cell(nsubjects,cfg.nPerm);


if ~exist(fullfile(outputDir,[cfg.outputName '.mat']),'file')
    for sub = 1:nsubjects
        
        fprintf('Processing subject %s (%d/%d) \n',cfg.subjects{sub},sub,length(cfg.subjects))
        
        % set random generator for repeatability
        rng(1,'twister')
        
        %% Get the data
        for p = 1:cfg.nPerm % shuffle labels
            
            load(fullfile(cfg.root,'Results',cfg.subjects{sub},'PrepData'),...
                'detection','stimulus','response');
            
            fprintf('\t Permutation %d of %d \n',p,cfg.nPerm+1)
            
            % select trials
            class1_idx = find(eval(cfg.conIdx{1}));
            class2_idx = find(eval(cfg.conIdx{2}));
            class_idx  = [class1_idx; class2_idx];
            nTrls      = length(class_idx);           
            
            % permute labels
                shuffleidx = randperm(nTrls);
            
            for per = 1:cfg.nPerm % random downsamples                
                
                load(fullfile(cfg.root,'Results',cfg.subjects{sub},'PrepData'));                
                
                % select trials                
                confidence = confidence(class_idx);
                data       = data(class_idx,:);
                detection  = detection(class_idx);
                response   = response(class_idx,:);
                run_idx    = run_idx(class_idx);
                stimulus   = stimulus(class_idx);
                                
                % shuffle
                confidence = confidence(shuffleidx);
                detection  = detection(shuffleidx);
                response   = response(shuffleidx,:);
                stimulus   = stimulus(shuffleidx);
                
                % determine labels
                labels     = eval(cfg.conIdx{2})+1;
                
                % Balance the confidence per run
                nRuns = length(unique(run_idx));
                
                trl_ind = [];
                for r = 1:nRuns
                    ind = find(run_idx == r);
                    lab_r = labels(ind);
                    ind1 = ind(lab_r==1);
                    ind2 = ind(lab_r==2);
                    for c = 1:6
                        cidx1 = ind1(confidence(ind1)==c);
                        cidx2 = ind2(confidence(ind2)==c);
                        if length(cidx1) > length(cidx2)
                            nDiff = length(cidx1)-length(cidx2);
                            cidx1 = cidx1(randperm(length(cidx1)));
                            cidx1(1:nDiff) = []; % randomly remove as many trls
                        elseif length(cidx2) > length(cidx1)
                            nDiff = length(cidx2)-length(cidx1);
                            cidx2 = cidx2(randperm(length(cidx2)));
                            cidx2(1:nDiff) = []; % randomly remove as many trls
                        end
                        
                        ctrls = [cidx1; cidx2];
                        
                        trl_ind = [trl_ind; ctrls]; clear ctrls cidx1 cidx2
                    end
                    
                end
                
                confidence = confidence(trl_ind);
                data       = data(trl_ind,:);
                detection  = detection(trl_ind);
                response   = response(trl_ind,:);
                run_idx    = run_idx(trl_ind);
                stimulus   = stimulus(trl_ind);
                
                Y          = eval(cfg.conIdx{1});
                nTrials    = length(Y);
                
                % zscore per run
                for r = 1:nRuns
                    data(run_idx==r,:) = zscore(data(run_idx==r,:),[],1);
                end
                
                %% Do decoding per ROI
                % decoding settings
                cfgD.gamma = cfg.gamma;
                
                % create folds - leave one run out
                folds = cell(nRuns,1);
                for r = 1:nRuns
                    folds{r} = find(run_idx==r);
                end
                nFolds = length(folds);
                
                % run over ROIs
                for r = 1:nROIs
                    
                    % mask the betas
                    x = data(:,idx{r}); % use mask indices to get data here
                    [~,b] = sort(abs(mean(x)));
                    if size(x,2) > cfg.nvox
                        x     = x(:,b(1:cfg.nvox)); % take n most active voxels
                    end
                    
                    % decoding
                    yhat = nan(length(Y),1);
                    for f = 1:nFolds
                        testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
                        labels = Y(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
                        
                        % train
                        decoder = train_LDA(cfgD,labels,trainX');
                        
                        % decode
                        yhat(testidx) = decode_LDA(cfgD,decoder,testX');
                    end
                    
                    % determine accuracy
                    permAcc(sub,r,per,p) = mean(yhat > 0 == Y);
                end
                
                
            end
            
        end
        
    end
    
    % save
    save(fullfile(outputDir,[cfg.outputName '.mat']),'permAcc');
else
    % load
    load(fullfile(outputDir,[cfg.outputName '.mat']),'permAcc');
    
end

%% Create bootstrap distributions
permAccAvg = squeeze(mean(permAcc,3)); % mean over different downsamples

bAccuracy_pm    = zeros(nROIs, cfg.nBtstrp);
for b =  1:cfg.nBtstrp
    
    if mod(b,1000)==0; fprintf('Bootstrapping %d/%d \n',b,cfg.nBtstrp); end
    
    tmp_pm = zeros(nsubjects,nROIs);
    for s = 1:nsubjects
        tmp_pm(s,:) = permAccAvg(s,:,randi(cfg.nPerm));%permAcc(s,:,randi(cfg.nPerm),randi(cfg.nPerm));
    end
    
    
    bAccuracy_pm(:,b) = mean(tmp_pm,1); clear tmp_pm
end

%% Plot the results

% get accurcy downsampled data
load(fullfile(cfg.root,'Results','GroupResults',cfg.outputDir,[cfg.outputName '.mat']),'accuracy');
load(fullfile(cfg.root,'Results','GroupResults',cfg.outputDir,[cfg.outputName '_randomdel.mat']),'permAcc_rd');

mAcc     = squeeze(mean(accuracy,3));
mAccRD   = squeeze(mean(permAcc_rd,3));

% just equal confidence versus null\
figure
for r = 1:length(roidx)
    subplot(2,2,r);
    histogram(bAccuracy_pm(roidx(r),:),'BinWidth',0.001,'FaceColor','r'); hold on;
    
    plot([mean(mAcc(:,roidx(r))) mean(mAcc(:,roidx(r)))],[0 1000],'y','LineWidth',2); hold on;
    plot([mean(mAccRD(:,roidx(r))) mean(mAccRD(:,roidx(r)))],[0 1000],'b','LineWidth',2); hold on;
    
    xlim([0.46 0.56]); ylim([0 1000])
    title(names{roidx(r)})
end
pValsAcc = sum(bAccuracy_pm(roidx,:)'>mean(mAcc(:,roidx),1))./cfg.nBtstrp
pValsAccRD = sum(bAccuracy_pm(roidx,:)'>mean(mAccRD(:,roidx),1))./cfg.nBtstrp

%% Plot them boxplots
map= makeColorMaps('dusk');
cs = map(round(linspace(20,length(map)-20,nROIs)),:);

figure;
for r = 1:length(roidx)
    subplot(1,4,r)
    dat = [mAcc(:,roidx(r)) mAccRD(:,roidx(r))];
    plot([1:2]+randn(length(cfg.subjects),1)*0.1,dat,'Color',cs(r,:),...
        'LineStyle','none','marker','*','Linewidth',2); hold on;
    b = boxplot(dat,'Colors','k','Symbol','r');
    hold on; set(b,{'linew'},{2});
    hold on; plot(xlim,[0.5 0.5],'k--','LineWidth',2);
    ylim([0.25 0.9])
    title(names{roidx(r)});
end

[~,pValsComp] = ttest(mAcc(:,roidx),mAccRD(:,roidx))
