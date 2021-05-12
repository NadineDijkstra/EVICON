function PrepareData(cfg)

nsubjects = length(cfg.subjects);

load(fullfile(cfg.outputDir,'general_mask'),'mask');

for sub = 1:nsubjects
    
    outDir = fullfile(cfg.outputDir,cfg.subjects{sub});
    if ~exist(outDir,'dir'); mkdir(outDir); end
    output = fullfile(outDir,'PrepData.mat');
    if ~exist(output,'file')
        fprintf('Preparing data for %s \n', cfg.subjects{sub});
    
    % get the trial info
    load(fullfile(cfg.dataDir,cfg.subjects{sub},'SPM.mat'));
    
    trial_indices = find(contains(SPM.xX.name,'detection') & ~contains(SPM.xX.name,'missed')); % trls in which we got a response
    
    nTrials = length(trial_indices);    
    detection  = nan(nTrials,1); 
    stimulus   = nan(nTrials,1);
    response   = nan(nTrials,2); % response and RT
    confidence = nan(nTrials,1);
    run_idx    = nan(nTrials,1);
    
    data       = zeros(nTrials,sum(mask(:)));
    for t = 1:nTrials
        
        [~,beta] = read_nii(fullfile(cfg.dataDir,cfg.subjects{sub},...
            sprintf('beta_%04d.nii',trial_indices(t))));
        data(t,:) = beta(mask);
        
        name = SPM.xX.name{trial_indices(t)};
        idx = strfind(name,'detection');
        detection(t) = str2double(name(idx+10));
        idx = strfind(name,'stimulus');
        stimulus(t) = str2double(name(idx+9));
        idx = strfind(name,'response');
        response(t,1) = str2double(name(idx+9));
        idx = strfind(name,'RT');
        response(t,2) = str2double(name(idx+3:idx+6));
        idx = strfind(name,'confidence');
        confidence(t) = str2double(name(idx+11));
        run_idx(t) = str2double(name(4));
    end
    
    % save the data
    save(output,'data','confidence','detection','stimulus','response','run_idx')
    clearvars -except sub cfg mask nsubjects
    else
        fprintf('Already prepared data for %s \n', cfg.subjects{sub});
    end 
    
end