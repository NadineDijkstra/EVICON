function bb_separate_clusters_cmdLine(cfg)
% bb_separate_clusters_cmdLine(cfg)
%
% This function distinguishes clusters within a statistics image and saves
% them in separate nii-files.
%
% asks user to manually put FSL in command line before continuing
%
% cfg.inputfile       = statisical image (e.g. t-map)
% cfg.threshold       = t-value
% cfg.numVox          = minimum num voxels per cluster (default = 0)

% load variables from cfg
get_vars_from_struct(cfg)


%%

[root, name, ext] = fileparts(inputfile);
cd(root);
name = strtok(name,'.');
inputFn  = [name '.nii.'];
outputFn = [name '_clusters.nii.gz'];
tableFn  = [name '_clusters.txt'];

fslCmd = sprintf('cluster -i %s --thresh=%f -o %s --mm > %s',inputFn,threshold,outputFn,tableFn);
fprintf('Cd to %s and run the following cmd and press any key when finished: \n \t %s \n',root,fslCmd);
pause;

fid = fopen(tableFn);
cellTable        = textscan(fid,'%s %s %s %s %s %s %s %s %s');
cellTable        = horzcat(cellTable{:});
cellTable(1:3,:) = [];
cellTable(:,7:9) = [];

numVoxels = str2double(cellTable(:,2));
clusterIdx = str2double(cellTable(:,1));
cIdx     = clusterIdx(numVoxels >= cfg.numVox);

% extract and save big enough clusters
gunzip(outputFn); delete(outputFn);
[V, Y] = read_nii(fullfile(root, [name '_clusters.nii']));

Y(~ismember(Y,cIdx)) = 0;

for i = 1:length(cIdx)
    Y(Y==cIdx(i)) = i;
end

write_nii(V,Y,fullfile(root, [name '_clusters.nii']));

% print big enough clusters
fprintf('\n Cluster index - numVoxels - rpVal - X - Y - Z \n')
disp(str2double(cellTable(ismember(clusterIdx,cIdx),:)))

