function DecodingSearchlight(cfg)
% function DecodingSearchlight(cfg)

% get the mask % searchlights indices
load(fullfile(cfg.root,cfg.mask),'mask','vind','mind'); 
V = read_nii('Results\general_mask.nii');

% loop over subjects
for sub = 1:length(cfg.subjects)
    
    fprintf('PROCESSING SUBJECT %s (%d/%d) \n',cfg.subjects{sub},sub,length(cfg.subjects))
    
    % set random generator for repeatability
    rng(1,'twister')
        
    
    %% Get the data
    outputDir = fullfile(cfg.root,'Results', cfg.subjects{sub});
    if ~exist(fullfile(outputDir,sprintf('Decoding_%s.mat',cfg.outputName)),'file')
    
    load(fullfile(outputDir,'PrepData'));
    
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
    idx = cell(nRuns,1); trl_ind = [];
    for r = 1:nRuns
        ind = find(run_idx == r);
        idx{r} = balance_trials(labels(ind,1),'downsample');
        
        trl_ind = [trl_ind; ind(cell2mat(idx{r}))];
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
    
    %% Do decoding per searchlight
    % decoding settings
    cfgD.gamma = cfg.gamma;
    ind        = find(mask);
    
    nSearchlights = length(vind);
    accuracy = zeros(V.dim);
    Yhat     = zeros(nTrials,nSearchlights);
    
    % create folds - leave one run out
    folds = cell(nRuns,1);
    for r = 1:nRuns
        folds{r} = find(run_idx==r);
    end
    nFolds = length(folds);
    
    % run over searchlights    
    for s = 1:nSearchlights
        
        if s >= (nSearchlights/10) && mod(s,round((nSearchlights/10))) == 0
            fprintf('\t Progress: %d percent of searchlights \n',round((s/nSearchlights)*100))
        end
        
        % mask the betas
        x = data(:,mind{s}); % use mask indices to get data here
        
        % decoding
        for f = 1:nFolds
            testidx = folds{f}; trainidx = setdiff(1:nTrials,testidx);
            labels = Y(trainidx); trainX = x(trainidx,:); testX = x(testidx,:);
            
            % train
            decoder = train_LDA(cfgD,labels,trainX');
            
            % decode
            Yhat(testidx,s) = decode_LDA(cfgD,decoder,testX');
        end
        
        % determine accuracy
        accuracy(ind(s)) = mean(Yhat(:,s) > 0' == Y);
        clear x
        
    end
    
    % write results
    write_nii(V, accuracy, fullfile(outputDir,...
        sprintf('Accuracy_%s.nii',cfg.outputName)))
    save(fullfile(outputDir,sprintf('Decoding_%s',cfg.outputName)),...
        'Yhat','Y','confidence','stimulus','response','data',...
        'detection','run_idx','-v7.3');    
    clear accuracy data confidence Yhat stimulus response detection run_idx
    else
        warning('Already processed this subject, skipping for now...')
    end
    
end
