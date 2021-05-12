function DecodingROIs(cfg)
% function DecodingROIs(cfg)

% outputDir
outputDir = fullfile(cfg.root,'Results', 'GroupResults',cfg.outputDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% get the ROI indices
load(fullfile(cfg.root,cfg.ROIs),'idx','names');

% loop over subjects
nsubjects  = length(cfg.subjects);
nROIs      = length(names);
accuracy   = zeros(nsubjects,nROIs);
permAcc    = zeros(nsubjects,nROIs,cfg.nPerm);
yhat       = cell(nsubjects,nROIs);
conf       = cell(nsubjects,1);
resp       = cell(nsubjects,1);
corrY      = cell(nsubjects,1);

if ~exist(fullfile(outputDir,[cfg.outputName '.mat']),'file')
    for sub = 1:nsubjects
        
        fprintf('Processing subject %s (%d/%d) \n',cfg.subjects{sub},sub,length(cfg.subjects))
        
        % set random generator for repeatability
        rng(1,'twister')
        
        %% Get the data
        load(fullfile(cfg.root,'Results',cfg.subjects{sub},'PrepData'));
        
        % select trials
        class1_idx = find(eval(cfg.conIdx{1}));
        class2_idx = find(eval(cfg.conIdx{2}));
        class_idx  = [class1_idx; class2_idx];
        
        confidence = confidence(class_idx);
        data       = data(class_idx,:);
        detection  = detection(class_idx);
        response   = response(class_idx,:);
        run_idx    = run_idx(class_idx);
        stimulus   = stimulus(class_idx);
        
        labels     = eval(cfg.conIdx{2})+1;
        
        % Balance the trials per run
        nRuns = length(unique(run_idx));
        idx2 = cell(nRuns,1); trl_ind = [];
        for r = 1:nRuns
            ind = find(run_idx == r);
            idx2{r} = balance_trials(labels(ind,1),'downsample');
            trl_ind = [trl_ind; ind(cell2mat(idx2{r}))];
        end
        
        confidence = confidence(trl_ind);
        data       = data(trl_ind,:);
        detection  = detection(trl_ind);
        response   = response(trl_ind,:);
        run_idx    = run_idx(trl_ind);
        stimulus   = stimulus(trl_ind);
        conf{sub}  = confidence;
        resp{sub}  = response;
        
        Y          = eval(cfg.conIdx{1});
        nTrials    = length(Y);
        corrY{sub} = Y;
        
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
            
            fprintf('\t Decoding ROI %d out of %d \n',r,nROIs)
            
            % mask the betas
            x = data(:,idx{r}); % use mask indices to get data here
            [~,b] = sort(abs(mean(x)));
            if size(x,2) > cfg.nvox
                x     = x(:,b(1:cfg.nvox)); % take n most active voxels
            end
            
            % decoding
            for f = 1:nFolds
                testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
                labels = Y(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
                
                % train
                decoder = train_LDA(cfgD,labels,trainX');
                
                % decode
                yhat{sub,r}(testidx) = decode_LDA(cfgD,decoder,testX');
            end           
            
            % determine accuracy
            accuracy(sub,r) = mean(yhat{sub,r}' > 0 == Y);            
            
            % permute
            for p = 1:cfg.nPerm
                
                %if mod(p,20)==0; fprintf('Permutation %d/%d \n',p,cfg.nPerm); end
                
                permY = zeros(length(Y),1);
                for f = 1:nFolds % permute within each run
                    labels = Y(folds{f});
                    permY(folds{f}) = labels(randperm(length(labels)));
                    clear labels
                end
                
                for f = 1:nFolds
                    testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
                    labels = permY(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
                    
                    % train
                    decoder = train_LDA(cfgD,labels,trainX');
                    
                    % decode
                    yhatP(testidx) = decode_LDA(cfgD,decoder,testX');
                end
                permAcc(sub,r,p) = mean(yhatP' > 0 == permY); clear yhatP           
            end
            clear x
            
        end
        
        
        
    end
    
    % save
    save(fullfile(outputDir,[cfg.outputName '.mat']),'accuracy','yhat',...
        'corrY','resp','conf','permAcc');

%% Do the second level bootstrapping

bAccuracy = zeros(nROIs, cfg.nBtstrp);
for b =  1:cfg.nBtstrp
    
   if mod(b,1000)==0; fprintf('Bootstrapping %d/%d \n',b,cfg.nBtstrp); end
   
   tmp = zeros(nsubjects,nROIs);
   for s = 1:nsubjects
       tmp(s,:) = permAcc(s,:,randi(cfg.nPerm)); % pick random perm
   end
   
   bAccuracy(:,b) = mean(tmp,1); clear tmp 
    
end

% append
save(fullfile(outputDir,[cfg.outputName '.mat']),'bAccuracy','-append');

else
    % load
    load(fullfile(outputDir,[cfg.outputName '.mat']),'accuracy','bAccuracy');
    
end

%% plot the results
figure;
map= makeColorMaps('dusk');
cs = map(round(linspace(20,length(map)-20,nROIs)),:);

pVals = sum(bAccuracy'>mean(accuracy,1))./cfg.nBtstrp

for r = 1:nROIs
    dat = squeeze(accuracy(:,r));
    plot(r+randn(length(cfg.subjects),1)*0.1,dat,'Color',cs(r,:),...
        'LineStyle','none','marker','*','Linewidth',2); hold on;    
    
    % significant?
    if pVals(r) < 0.05; hold on; plot(r,0.8,'k*'); end
end
b = boxplot(accuracy,'Colors','k','Symbol','r');
hold on; set(b,{'linew'},{2});
hold on; plot(xlim,[0.5 0.5],'k--','LineWidth',2);
ylim([.25 0.9])
set(gca,'XTick',1:nROIs); set(gca,'XTickLabels',names);

ylim([.25 0.9])

