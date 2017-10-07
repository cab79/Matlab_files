function Combine_clusters_source(S)

%% Requires input structure S containing fields as follows. See Cluster_processing script for examples.
%-------------------------------------------------------------
% S.spmstats_path: directory in which SPM analysis is saved 
% S.spm_dir: specific folder containing the SPM stats for this analysis
% S.contrasts: contrast name cell array - must match that in Matlabbatch (i.e. from design-batch script). Leave empty to proccess ALL contrasts in Matlabbatch
% S.clustab{tf(cont)}: stats headers to save in a table 
% S.batch: name of batch .mat file saved from design_batch.m and within same folder as SPM.mat

%OPTIONAL SETTINGS (advise not to change):
%Num=16;    %- number of maxima per cluster [3]
%Dis=8;    %- min distance among clusters {mm} [8]
%Str='';    %- header string

%% RUN
if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

% load all c*_spm.nii files from all contrasts of interest and store in a
% temporary cell array
niicell = {};
n=0;
for sp = 1:length(S.spm_paths)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp});

    if isempty(S.contrasts)
        % load SPM design and list contrasts
        load(fullfile(S.spm_path,S.batch));
        contrasts = {};
        tf=[];
        for fc = 1:length(matlabbatch{3}.spm.stats.con.consess)
            try
                contrasts{fc} = matlabbatch{3}.spm.stats.con.consess{1,fc}.fcon.name;
                tf=[tf 1];
            catch
                contrasts{fc} = matlabbatch{3}.spm.stats.con.consess{1,fc}.tcon.name;
                tf=[tf 2];
            end
        end
    else
        contrasts=S.contrasts;
        tf=S.tf;
    end

    for cont = 1:length(contrasts)
        contrast = contrasts{cont};
        conname = strrep(contrast,' ','');
        conname = strrep(conname,'*','_');
        clus_path = fullfile(S.spm_path,[conname '_clusters']); 
        
        cfnames = dir(fullfile(clus_path,'*_spm.nii'));
        for i = 1:length(cfnames)
            % load each nii img array into a cell
            n = n+1;
            nii = load_nii(fullfile(S.clus_path,cfname(i).name));
            niicell{n} =  nii.img;
        end

    end
end
        
% for each cell, cycle through to identify unique voxel value in new image
comnii = zeros(size(niicell{1})); % new image
for n = 1:length(niicell)
    % current max value in comnii
    maxval = max(comnii);
    % adds new values that are unique to comnii
    comnii = comnii+(maxval+1)*reshape((niicell{n}>1),size(comnii));
end
        
% save the new image
imgpath = fullfile(S.spmstats_path,'Combined_clusters');
if ~isdir(imgpath)
    mkdir(imgpath)
end

save_nii(comnii,fullfile(imgpath,'comb_clus'));

