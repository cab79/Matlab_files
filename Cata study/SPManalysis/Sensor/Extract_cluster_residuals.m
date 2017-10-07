function Extract_cluster_residuals(S)
%% for extracting cluster residuals, everaged over time and space, for normality tests.
dbstop if error
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
voi_pref = 'VOI_'; %VOI filename(s) prefix

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
    imglist = dir(fullfile(S.spm_path,['Res_*.nii'])); % Subject-condition residual image list
    S.imglist = {imglist(:).name};

    S.clus_path={};
    %if ~isfield(S,'clus_path')
        alldir=dir(fullfile(S.spm_path,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(S.spm_path,alldir(c).name);
        end
    %end

    for cldir = 1:length(S.clus_path)

        cimages = dir(fullfile(S.clus_path{cldir},gclusname));

        %load VOI filenames
        %vfiles = dir(fullfile(S.clus_path{cldir},[voi_pref '*.xlsx']));

        for c = 1:length(cimages)
            [~,cname,~] = fileparts(cimages(c).name);

            % load VOI excel file (as template for new residuals data)
            voiname = strsplit(cname,'_');
            vname = [voi_pref voiname{1} '.xlsx'];
            [~,~,res_temp] = xlsread(fullfile(S.clus_path{cldir},vname));
            dat_col = find(strcmp(res_temp(1,:),'Data'));

            % load cluster image
            Cnii = load_nii(fullfile(S.clus_path{cldir},cimages(c).name));

            % for each data image file,
            res=[];
            for i = 1:length(S.imglist)
                disp(['contrast ' num2str(cldir) '/' num2str(length(S.clus_path)) ', cluster ' num2str(c) '/' num2str(length(cimages)) ', image ' num2str(i) '/' num2str(length(S.imglist))]); % display progress
                iname = S.imglist{i};
                nii = load_nii(fullfile(S.spm_path,iname));
                im = nii.img.*Cnii.img;
                res(i,1)=nanmean(im(:),1);
            end

            res_temp(2:end,dat_col) = num2cell(res);

            % save 
            sname = strrep(vname,'VOI','Res');
            xlswrite(fullfile(S.clus_path{cldir},sname),res_temp);

        end

    end
end