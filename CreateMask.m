function CreateMask(cfg)

v2struct(cfg); % unpack 

[V,gm_mask] = read_nii(gm_mask);
m = [];
for sub = 1:length(subjects)
    [~,map] = read_nii(fullfile('Data',subjects{sub},'beta_0001.nii'));
    m = cat(4,m,map); clear map
end
nan_mask = any(isnan(m),4);    
mask = ~nan_mask & gm_mask > 0.1;

% write mask as nifti
write_nii(V,double(mask),'Results\general_mask.nii');

% get searchlight indices based on mask
[vind,mind] = searchlightIndices(mask,searchlight);

save('Results\general_mask.mat','mask','vind','mind');
