function DecodingROIsDistrBalance(cfg)
% function DecodingROIsDistrBalance(cfg)

% outputDir
outputDir = fullfile(cfg.root,'Results','GroupResults',cfg.outputDir);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% get the ROI indices
load(fullfile(cfg.root,cfg.ROIs),'idx','names');

% loop over subjects
nsubjects  = length(cfg.subjects);
nROIs      = length(names);
accuracy   = zeros(nsubjects,nROIs,cfg.nPerm);
yhat       = cell(nsubjects,nROIs,cfg.nPerm,cfg.nPerm);
conf       = cell(nsubjects,cfg.nPerm);
resp       = cell(nsubjects,cfg.nPerm);
corrY      = cell(nsubjects,cfg.nPerm);
ridx       = cell(nsubjects,cfg.nPerm);

if ~exist(fullfile(outputDir,[cfg.outputName '.mat']),'file')
    for sub = 1:nsubjects
        
        fprintf('Processing subject %s (%d/%d) \n',cfg.subjects{sub},sub,length(cfg.subjects))
        
        % set random generator for repeatability
        rng(1,'twister')
        
        %% Get the data
        for per = 1:cfg.nPerm
            
            
            load(fullfile(cfg.root,'Results',cfg.subjects{sub},'PrepData'));        
            fprintf('\t Permutation %d of %d \n',per,cfg.nPerm)
        
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
            conf{sub,per}  = confidence;
            resp{sub,per}  = response;
            ridx{sub,per}  = run_idx;
            
            Y          = eval(cfg.conIdx{1});
            nTrials    = length(Y);
            corrY{sub,per} = Y;
            
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
                
                %fprintf('\t Decoding ROI %d out of %d \n',r,nROIs)
                
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
                    yhat{sub,r,per}(testidx) = decode_LDA(cfgD,decoder,testX');
                end
                
                % determine accuracy
                accuracy(sub,r,per) = mean(yhat{sub,r,per}' > 0 == Y);            
                
            end    
            
            
        end
        
    end
    
    % save
    save(fullfile(outputDir,[cfg.outputName '.mat']),'accuracy','yhat',...
        'corrY','resp','conf','ridx');
else
    % load
    load(fullfile(outputDir,[cfg.outputName '.mat']),'accuracy','corrY','conf','permAcc','ridx');
    
end

%% Now randomly delete as many trls from each class to see if it's a power thing
% Generate null-distribution based on this

% loop over subjects
permAcc_rd    = zeros(nsubjects,nROIs,cfg.nPerm);
if ~exist(fullfile(outputDir,[cfg.outputName '_randomdel.mat']),'file')
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
        
        Confidence = confidence(class_idx);
        Data       = data(class_idx,:);
        detection  = detection(class_idx);
        Detection  = detection;
        response   = response(class_idx,:);
        Response   = response;
        Run_idx    = run_idx(class_idx);
        stimulus   = stimulus(class_idx);
        Stimulus   = stimulus;
        
        Labels     = eval(cfg.conIdx{2})+1;
        
        % Loop over permutations
        % Downsample the trials per run
        nRuns = length(unique(run_idx));
        for p = 1:cfg.nPerm
            fprintf('\t Calculatin pemrutation %d out of %d \n',p,cfg.nPerm)
            
            trl_ind = [];
            for r = 1:nRuns
                ind   = find(Run_idx == r);
                nc1   = sum(corrY{sub}(ridx{sub}==r)==1);
                nc2   = sum(corrY{sub}(ridx{sub}==r)==0);
                
                ind1  = ind(Labels(ind)==1);
                ind1  = ind1(randperm(length(ind1)));
                ind1  = ind1(1:nc1); % randomly select same number of trials
                ind2  = ind(Labels(ind)==2);
                ind2  = ind2(randperm(length(ind2)));
                ind2  = ind2(1:nc2);
                
                ind     = [ind1; ind2];
                trl_ind = [trl_ind; ind];
            end
            
            confidence = Confidence(trl_ind);
            data       = Data(trl_ind,:);
            detection  = Detection(trl_ind);
            response   = Response(trl_ind,:);
            run_idx    = Run_idx(trl_ind);
            stimulus   = Stimulus(trl_ind);
            conf_rd{sub,p} = confidence;
            
            Y          = eval(cfg.conIdx{1});
            nTrials    = length(Y);
            corrY_rd{sub,p} = Y;
            
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
                for f = 1:nFolds
                    testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
                    labels = Y(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
                    
                    % train
                    decoder = train_LDA(cfgD,labels,trainX');
                    
                    % decode
                    Yhat(testidx) = decode_LDA(cfgD,decoder,testX');
                end
                
                % determine accuracy
                permAcc_rd(sub,r,p) = mean(Yhat' > 0 == Y); clear Yhat
                
                clear x Yhat
            end
        end
    end
    
    % save
    save(fullfile(outputDir,[cfg.outputName '_randomdel.mat']),'conf_rd','corrY_rd','permAcc_rd');
else
    % load
    load(fullfile(outputDir,[cfg.outputName '_randomdel.mat']),'conf_rd','corrY_rd','permAcc_rd');
end

%% Plot the confidence distributions
cf_ed = zeros(6,2,2);
cf_rd = zeros(6,2,2);
for c = 1:6
    tmp_ed = zeros(nsubjects,2);
    tmp_rd = zeros(nsubjects,2);
    for s = 1:nsubjects
        tmp_ed(s,1) = sum(conf{s}(corrY{s}==1)==c)/sum(corrY{s}==1);
        tmp_ed(s,2) = sum(conf{s}(corrY{s}==0)==c)/sum(corrY{s}==0);
        
        tmp_rd_p = zeros(cfg.nPerm,2);
        for p = 1:cfg.nPerm
            tmp_rd_p(p,1) = sum(conf_rd{s,p}(corrY_rd{s,p}==1)==c)/sum(corrY_rd{s,p}==1);
            tmp_rd_p(p,2) = sum(conf_rd{s,p}(corrY_rd{s,p}==0)==c)/sum(corrY_rd{s,p}==0);
        end
        tmp_rd(s,1) = mean(tmp_rd_p(:,1),1);
        tmp_rd(s,2) = mean(tmp_rd_p(:,2),1);
    end
    
    for cl = 1:2
        cf_ed(c,cl,1) = mean(tmp_ed(:,cl));
        cf_ed(c,cl,2) = std(tmp_ed(:,cl))./sqrt(nsubjects);
        
        cf_rd(c,cl,1) = mean(tmp_rd(:,cl));
        cf_rd(c,cl,2) = std(tmp_rd(:,cl))./sqrt(nsubjects);
    end
end

figure;
subplot(2,1,1);
barwitherr(squeeze(cf_ed(:,:,2)),squeeze(cf_ed(:,:,1)));
xlabel('Confidence'); ylabel('Proportion'); legend('Class 1','Class 2');
title('Equal confidence');
subplot(2,1,2);
barwitherr(squeeze(cf_rd(:,:,2)),squeeze(cf_rd(:,:,1)));
xlabel('Confidence'); ylabel('Proportion'); legend('Class 1','Class 2');
title('Random downsample');
