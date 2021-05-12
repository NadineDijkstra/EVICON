% Set paths etc.
clear all;
restoredefaultpath;

root = 'D:\EVICON\';
cd(root)

addpath('Analyses');
addpath('Analyses\Utilities');
addpath('Analyses\spm12'); % parth to your SPM directory

% Subjects
subjects = {'sub-MP02312', 'sub-MP02319', 'sub-MP02320', 'sub-MP02324', 'sub-MP02325',...
    'sub-MP02327', 'sub-MP02328', 'sub-MP02332', 'sub-MP02333', 'sub-MP02336',...
    'sub-MP02337', 'sub-MP02338', 'sub-MP02339', 'sub-MP02341', 'sub-MP02348',...
    'sub-MP02350', 'sub-MP02351', 'sub-MP02353', 'sub-MP02360', 'sub-MP02361',...
    'sub-MP02363', 'sub-MP02364', 'sub-MP02376', 'sub-MP02377', 'sub-MP02378',...
    'sub-MP02399', 'sub-MP02407', 'sub-MP02417', 'sub-MP02425', 'sub-MP02440',...
    'sub-MP02500', 'sub-MP02504', 'sub-MP02522', 'sub-MP02546', 'sub-MP02552'};
nsubjects = length(subjects);

%% 1. Prepare data
% 1.1 Create mask based on no nans in all subs and grey matter mask
cfg = [];
cfg.root        = root;
cfg.subjects    = subjects; 
cfg.gm_mask     = 'Results\rgm_mask.nii';
cfg.searchlight = 4; 

CreateMask(cfg)

% 1.2 Extract trial information from data and select voxels within mask 
cfg = [];
cfg.subjects = subjects;
cfg.dataDir  = fullfile(root,'Data');
cfg.outputDir= fullfile(root,'Results');
PrepareData(cfg);

%% 2. Decode presence and identity searchlight

cfg = [];
cfg.root       = root;
cfg.subjects   = subjects;
cfg.mask       = '/Results/general_mask.mat';
cfg.gamma      = 0.1;

% 2.1 Discrimination
cfg.outputName = 'Discrimination_Corr'; 
cfg.conIdx{1} = "detection==0 & stimulus==1 & response(:,1)==1"; % class 1
cfg.conIdx{2} = "detection==0 & stimulus==3 & response(:,1)==3"; % class 2
DecodingSearchlight(cfg)

% 2.2 Detection
cfg.outputName = 'Detection_Corr'; 
cfg.conIdx{1} = "detection==1 & stimulus==0 & response(:,1)==0"; % class 1
cfg.conIdx{2} = "detection==1 & ismember(stimulus,[1,3]) & ismember(response(:,1),[1,3])"; % class 2
DecodingSearchlight(cfg)


%% 3. Decoding whole-brain stats 

% get permutation maps 
cfg = [];
cfg.root       = root;
cfg.subjects   = subjects;
cfg.mask       = '/Results/general_mask.mat';
cfg.nPerm      = 25; % per subject - combine to create group null-distr
cfg.gamma      = 0.1;

% 3.1 Discrimination
cfg.outputName = 'Discrimination_Corr'; 
cfg.conIdx{1} = "detection==0 & stimulus==1 & response(:,1)==1"; % class 1
cfg.conIdx{2} = "detection==0 & stimulus==3 & response(:,1)==3"; % class 2
DecodingSearchlightPermutation(cfg)

% 3.2 Detection
cfg.outputName = 'Detection_Corr'; 
cfg.conIdx{1} = "detection==1 & stimulus==0 & response(:,1)==0"; % class 1
cfg.conIdx{2} = "detection==1 & ismember(stimulus,[1,3]) & ismember(response(:,1),[1,3])"; % class 2
DecodingSearchlightPermutation(cfg)

% 3.3 Create permutation group distribution
cfg = [];
cfg.root       = root;
cfg.subjects   = subjects;
cfg.mask       = '/Results/general_mask.nii';
cfg.contrast   = 'Discrimination_Corr';

GroupDecodingPermutation(cfg)

% 3.4 FDR correction for multiple comparisons
cfg.inputfile           = fullfile(root,'Results','GroupResults',['pval_' cfg.contrast '.nii']);
cfg.qvalue            = 0.05; % FDR threshold (default = 0.05)
cfg.mask              = fullfile(root,'Results','general_mask.nii');

[V,map] = read_nii(cfg.inputfile); [~,mask] = read_nii(cfg.mask);
[~,fdr_threshold] = fdr_bh(map(mask>0),cfg.qvalue);    

% 3.5 get cluster labels
cfg.inputfile     = fullfile(root,'Results','GroupResults',['rpval_' cfg.contrast '.nii']);
cfg.threshold     = 1-fdr_threshold; % above 0
cfg.numVox        = 50;

bb_separate_clusters_cmdLine(cfg);

% 3.6 mask accuracy image significance
contrast = 'Discrimination_Corr';

[V,acc] = read_nii(fullfile(root,'Results','GroupResults',['Accuracy_' contrast '.nii']));
[V,clusters] = read_nii(fullfile(root,'Results','GroupResults',['rpval_' contrast '_clusters.nii']));
mAcc =  acc; mAcc(~clusters>0) = 0;
V.dt = [4 0]; V.pinfo = [0 0.5 352]';
write_nii(V,mAcc,fullfile(root,'Results','GroupResults',['sig_acc_' contrast  '.nii']));


%% 4. Decode presence and identity from ROIs

cfg = [];
cfg.root       = root;
cfg.subjects   = subjects;
cfg.ROIs       = '/Results/ROI_mask_indices.mat';
cfg.gamma      = 0.1;
cfg.nvox       = 200;
cfg.nPerm      = 25; 
cfg.nBtstrp    = 10000; 
cfg.outputDir  = 'ROI_decoding';

% 4.1 Discrimination
cfg.outputName = 'Discrimination_Corr'; 
cfg.conIdx{1} = "detection==0 & stimulus==1 & response(:,1)==1"; % class 1
cfg.conIdx{2} = "detection==0 & stimulus==3 & response(:,1)==3"; % class 2
DecodingROIs(cfg)

% 4.2 Detection
cfg.outputName = 'Detection_Corr'; 
cfg.conIdx{1} = "detection==1 & stimulus==0 & response(:,1)==0"; % class 1
cfg.conIdx{2} = "detection==1 & ismember(stimulus,[1,3]) & ismember(response(:,1),[1,3])"; % class 2
DecodingROIs(cfg)

%% 5. Downsample confidence decoding
cfg = [];
cfg.root       = root;
cfg.subjects   = subjects;
cfg.ROIs       = '/Results/ROI_mask_indices.mat';
cfg.gamma      = 0.1;
cfg.nvox       = 200;
cfg.outputDir  = 'ROI_decoding/Balanced/DistributionEqual';

% 5.1 Discrimination
cfg.outputName = 'Discrimination_Corr'; 
cfg.conIdx{1} = "detection==0 & stimulus==1 & response(:,1)==1"; % class 1
cfg.conIdx{2} = "detection==0 & stimulus==3 & response(:,1)==3"; % class 2
DecodingROIsDistrBalance(cfg)

cfg.nPerm      = 25; % create random permutations for stats
cfg.nBtstrp    = 10000; 
DecodingROIsDistrBalancePermutation(cfg)

% 5.2 Detection
cfg.outputName = 'Detection_Corr'; 
cfg.conIdx{1} = "detection==1 & stimulus==0 & response(:,1)==0"; % class 1
cfg.conIdx{2} = "detection==1 & ismember(stimulus,[1,3]) & ismember(response(:,1),[1,3])"; % class 2
DecodingROIsDistrBalance(cfg)

cfg.nPerm      = 25; % create random permutations for stats
cfg.nBtstrp    = 10000; 
DecodingROIsDistrBalancePermutation(cfg)

