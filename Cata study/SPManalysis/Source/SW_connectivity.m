function SW_connectivity(S)

addpath(genpath('C:\Data\Matlab\osl-core'))
addpath('C:\Data\Matlab\MEG-ROI-nets\ROInets')

% PREREQUISITES: 
% - An SPM.mat file containing a contrast on source-space ERP data 
% - Associated cluster data created from "Extract_clusters_source.m"

%% Requires input structure S containing fields as follows. See Cluster_processing script for examples.
%-------------------------------------------------------------
% S.data_path: root directory in which subject-specific folders are located
% S.spmstats_path: directory in which SPM analysis is saved 
% S.spm_dir: specific folder containing the SPM stats for this analysis
% S.contrasts: contrast name cell array - must match that in Matlabbatch (i.e. from design-batch script). Leave empty to proccess ALL contrasts in Matlabbatch
% S.clus_path: cell array of full paths to clusters to be analysed
%generic cluster data name
gdataname = 'cluster_data*.mat';

%% run

if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

% For each analysis (time window)
%last_spm_path='';
last_tw=0;
for sp = 1:size(S.spm_paths,1)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    
    % run analysis on combined clusters?
    if strcmp(S.gclusname,'comb_clus.nii')
        anapath = fullfile(S.spmstats_path,['Timewin_' num2str(S.spm_paths{sp,2})]);
        S.contrasts = {'Combined'};
        % continue if this timewindow was just analysis
        if last_tw==S.spm_paths{sp,2}
            continue
        end
        last_tw = S.spm_paths{sp,2};
    else
        anapath = S.spm_path;
    end
    
    %S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    %if strcmp(last_spm_path,S.spm_path)
    %    continue
    %end
    %last_spm_path = S.spm_path;

    load(fullfile(S.spm_path,'SPM.mat'));
    S.Fm = SPM.xX.I; % factor matrix
    S.imglist = SPM.xY.P; % Subject-condition image list
    
    S.clus_path={};
    if isempty(S.contrasts)
        alldir=dir(fullfile(anapath,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(anapath,alldir(c).name);
        end
    else
        for c = 1:length(S.contrasts)
            S.clus_path{c} = fullfile(anapath,[S.contrasts{c} '_clusters']);
        end
    end

    % For each contrast
    for cldir = 1:length(S.clus_path)
        S.cldir=cldir;

        cdata = dir(fullfile(S.clus_path{cldir},gdataname));

        % For each cluster
        for c = 1:length(cdata)
            disp(['*** Cluster = ' num2str(c) ' / ' num2str(length(cdata))])

            [~,cname,~] = fileparts(cdata(c).name);

            % load
            Sc = load(fullfile(S.clus_path{cldir},cdata(c).name));
            
            % update S with data from Sc
            f = fieldnames(Sc.S);
            for i = 1:length(f)
                if ~isfield(S,f{i})
                    S.(f{i}) = Sc.S.(f{i})
                end
            end
            
            % Get cluster names
            fnames = fieldnames(S.wf);
            
            % only analyse first one for now
            wfall=S.wf.(fnames{1}).wf;
            
            [tol_ord,S] = find_ranks(S,wfall,S.gclusname,S.clus_path{cldir});
            save(fullfile(S.clus_path{cldir},'tol_ord.mat'),'tol_ord');
            
            tryagain=1;
            while tryagain
                [orth_wf,ucol,tol_i] = find_ortho_symm_v3(S,wfall,tol_ord);

                % save 
                S.wf.(fnames{1}).orth_wf = orth_wf;
                S.wf.(fnames{1}).ucol = ucol;
                save(fullfile(S.clus_path{cldir},[cname '.mat']),'S');

                %if length(unique(ucol))~=length(ucol)
                    create_reduced_image(S,S.gclusname,S.clus_path{cldir},1,[]);
                %end

                % select data
                time = S.wf.comb_clus.time;
                subs = unique(S.Fm(:,1+S.subrow));

                % correlation analysis
                try
                    mats = source_correlation_analysis(S,subs,time,orth_wf);
                    tryagain=0;
                catch ME
                    if strcmp(ME.identifier,'dp_glasso:NotPosDefCovInput')
                        tol_ord = tol_ord(tol_i+1:end); 
                    end
                end
            end
            
        end
        
        % reformat results - correlationMats is a cell array of frequency bands
        S.nFreqBands = length(S.frequencyBands);
        S.nSessions = length(subs);  
        correlationMats = reformat_results(mats, S);

        %% Subject-level analysis to average over sessions in a fixed-effects manner
        % will be same output as first-level if there is one session per subject
        correlationMats = do_subject_level_glm(correlationMats, S);

        %% Group-level analysis
        % Find whole group means
        if strcmpi(S.paradigm, 'rest'),
            correlationMats = do_group_level_statistics(correlationMats, S);
        end%if

        % Perform group-level GLM
        if ~isempty(S.GroupLevel),
            %correlationMats = select_regions(correlationMats,[1]);
            correlationMats = do_group_level_glm_CAB(correlationMats, S);
        end%if
        
        % save 
        S.wf.(fnames{1}).correlationMats = correlationMats;
        save(fullfile(S.clus_path{cldir},[cname '.mat']),'S');

    end
end


end

function CM = select_regions(CM,rg)
    if ~isempty(rg)
        levnames = {'firstLevel','subjectLevel'};
        for ln = 1:length(levnames)
            for lnn = 1:length(CM{1}.(levnames{ln}))
                fnames = fieldnames(CM{1}.(levnames{ln})(lnn).cope);
                for f = 1:length(fnames)
                    dat = CM{1}.(levnames{ln})(lnn).cope.(fnames{f});
                    if isnumeric(dat) && length(size(dat))==3
                        CM{1}.(levnames{ln})(lnn).cope.(fnames{f}) = dat(rg,:,:);
                    end
                end
            end
        end
    end
end