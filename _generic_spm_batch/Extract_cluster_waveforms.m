function Extract_cluster_waveforms(S)
%% for extracting cluster activity over time (averaged over space) for plotting

% PREREQUISITES: 
% - An SPM.mat file containing a contrast on sensor-space ERP data 
% - Associated cluster data created from "Extract_clusters.m"

%% Requires input structure S containing fields as follows. See Cluster_processing script for examples.
%-------------------------------------------------------------
% S.data_path: root directory in which subject-specific folders are located
% S.spmstats_path: directory in which SPM analysis is saved 
% S.spm_dir: specific folder containing the SPM stats for this analysis
% S.contrasts: contrast name cell array - must match that in Matlabbatch (i.e. from design-batch script). Leave empty to proccess ALL contrasts in Matlabbatch
% S.clus_path: cell array of full paths to clusters to be analysed
%generic cluster image name
gclusname = 'c*spm.nii';

%% run

if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

for sp = 1:length(S.spm_paths)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp});

    load(fullfile(S.spm_path,'SPM.mat'));
    S.Fm = SPM.xX.I; % factor matrix
    S.imglist = SPM.xY.P; % Subject-condition image list

    S.clus_path={};
    %if ~isfield(S,'clus_path')
        alldir=dir(fullfile(S.spm_path,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(S.spm_path,alldir(c).name);
        end
    %end

    for cldir = 1:length(S.clus_path)

        cimages = dir(fullfile(S.clus_path{cldir},gclusname));

        S.wf=struct; % clears this variable for the next iteration
        for c = 1:length(cimages)

            [~,cname,~] = fileparts(cimages(c).name);

            % create empty cell array for waveform data
            S.wf.(cname) = cell(length(S.imglist),1); 

            % load cluster image
            Cnii = load_nii(fullfile(S.clus_path{cldir},cimages(c).name));

            % create maximum intensity projection (MIP) over time
            mip = max(Cnii.img,[],3);

            % normalise to 1
            mip = mip/max(mip(:));

            % for each data image file,
            for i = 1:length(S.imglist)
                disp(['contrast ' num2str(cldir) '/' num2str(length(S.clus_path)) ', cluster ' num2str(c) '/' num2str(length(cimages)) ', image ' num2str(i) '/' num2str(length(S.imglist))]); % display progress
                iname = S.imglist{i};
                iname = iname(1:end-2);
                nii = load_nii(iname);

                % at each time point, multiply 2D data image with cluster MIP and average over voxels to
                % create a cluster waveform
                Ndp=size(nii.img,3);
                S.wf.(cname){i} = NaN(Ndp,1);
                for dp = 1:Ndp
                    im = nii.img(:,:,dp).*mip;
                    S.wf.(cname){i}(dp)=nanmean(im(:),1);
                end
            end
        end

        % save 
        save(fullfile(S.clus_path{cldir},'cluster_data.mat'),'S');
    end
end