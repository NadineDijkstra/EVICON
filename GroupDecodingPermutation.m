function GroupDecodingPermutation(cfg)

% get mask
[V,mask] = read_nii(fullfile(cfg.root,cfg.mask));
nVox     = sum(mask(:)>0);

% load accuracy map
acc   = zeros(length(cfg.subjects),nVox);
for s = 1:length(cfg.subjects)
    [~,tmp] = read_nii(fullfile(cfg.root,'Results',cfg.subjects{s},...
        ['Accuracy_' cfg.contrast '.nii']));
    acc(s,:) = tmp(mask>0); clear tmp
end

mAcc = squeeze(mean(acc,1)); clear acc; % mean accruacy

% load permutations
nPerm = 25;
permutations = zeros(length(cfg.subjects),nPerm,nVox);
for s = 1:length(cfg.subjects)
    load(fullfile(cfg.root,'Results',cfg.subjects{s},...
        ['Perm_' cfg.contrast '.mat']),'accuracy')
    permutations(s,:,:) = accuracy; clear accuracy
end
    
% create group distribution
nBoot = 10000;
bAcc = zeros(nBoot,nVox);
for b = 1:nBoot
    
    if mod(b,100) == 0
        fprintf('Bootstrapping: %d / %d \n', b, nBoot);
    end
    
    % pick random map per sub
    tmp = zeros(length(cfg.subjects),nVox);
    for s = 1:length(cfg.subjects)
        tmp(s,:) = squeeze(permutations(s,randi(nPerm),:));
    end
    bAcc(b,:) = mean(tmp,1); % average
    
end

% calculate p_vals
pVal = nan(size(mask));
pVal(mask>0) = sum(bAcc>mAcc)/nBoot;

% save as image
write_nii(V,pVal,fullfile(cfg.root,'Results',...
    'GroupResults',['pval_' cfg.contrast '.nii']))
write_nii(V,1-pVal,fullfile(cfg.root,'Results',...
    'GroupResults',['rpval_' cfg.contrast '.nii']))