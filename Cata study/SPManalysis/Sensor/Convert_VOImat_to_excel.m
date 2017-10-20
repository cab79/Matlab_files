function Convert_VOImat_to_excel(S)
%% loads VOIs saved from SPM and convert to long-format excel file

% PREREQUISITES: 
% - An SPM.mat file containing a factor matrix for the analysis
% - Associated cluster data created from "Extract_clusters.m"
% - Associated cluster waveform data created from "Extract_cluster_waveforms.m"

% S.subinfo: name of subject information file in SPM directory (produced by design_batch)
% S.factlev: factor and level information. 
    %1: factor name in design matrix, 2: output factor name 3: factor levels. 
    % Factors can be in order of desired output, not necessarily of
    % input into SPM design matrix. Levels must be in the same order as SPM
    % design, but characters don't need to match anything.
% S.subrow = 4; % row of above factlev containing the subject factor
% S.batch: name of batch .mat file saved from design_batch.m and within same folder as SPM.mat

voi_pref = 'VOI_c'; %VOI filename(s) prefix

%% run

if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

for sp = 1:size(S.spm_paths,1)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});

    if ~isfield(S,'Fm');
        load(fullfile(S.spm_path,'SPM.mat'));
        S.Fm = SPM.xX.I; % factor matrix
    end

    % load SPM design and list factors
    load(fullfile(S.spm_path,S.batch));
    S.fact = {matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fac(:).name};

    %load subject information
    load(fullfile(S.spm_path,S.subinfo));
    S.factlev{S.subrow,3} = strcat('''',subID);

    % identify file indices relating to all factors
    fact_col = [];
    for c = 1:length(S.factlev) % using a loop ensures that fact_col is in correct order for condition labels applied later
        fact_col(c) = find(ismember(S.fact,S.factlev{c,1}));
    end
    fact_ind = S.Fm(:,1+fact_col);
    fact_lev = {};
    for c = 1:length(S.factlev)
        factname = S.factlev{c,2};
        fact_lev{1,c} = factname{:};
        fact_lev(2:1+length(fact_ind),c) = S.factlev{c,3}(fact_ind(:,c))';
    end

    %if ~isfield(S,'clus_path')
        alldir=dir(fullfile(S.spm_path,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(S.spm_path,alldir(c).name);
        end
    %end

    for cldir = 1:length(S.clus_path)

        %load VOI filenames
        vfiles = dir(fullfile(S.clus_path{cldir},[voi_pref '*.mat']));

        for v = 1:length(vfiles)
            disp(['contrast ' num2str(cldir) '/' num2str(length(S.clus_path)) ', cluster ' num2str(v) '/' num2str(length(vfiles))]); % display progress
            vname = vfiles(v).name;
            load(fullfile(S.clus_path{cldir},vname))
            Yxl = fact_lev;
            dat_col = size(fact_lev,2)+1;
            Yxl{1,dat_col} = 'Data';
            Yxl(2:1+length(Y),dat_col) = num2cell(Y);
            [~,sname,~] = fileparts(vname);
            xlswrite(fullfile(S.clus_path{cldir},[sname '.xlsx']),Yxl);
        end

        % save 
        %save(fullfile(S.clus_path{cldir},'cluster_data.mat'),'S');
    end
end