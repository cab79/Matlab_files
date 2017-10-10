function SW_connectivity(S)
%% for extracting cluster activity over time (averaged over space) for plotting

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
%generic cluster image name
gclusname = 'cluster_data*.mat';

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
last_spm_path='';
for sp = 1:size(S.spm_paths,1)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    if strcmp(last_spm_path,S.spm_path)
        continue
    end
    last_spm_path = S.spm_path;

    load(fullfile(S.spm_path,'SPM.mat'));
    S.Fm = SPM.xX.I; % factor matrix
    S.imglist = SPM.xY.P; % Subject-condition image list
    
    S.clus_path={};
    if isempty(S.contrasts)
        alldir=dir(fullfile(S.spm_path,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(S.spm_path,alldir(c).name);
        end
    else
        for c = 1:length(S.contrasts)
            S.clus_path{c} = fullfile(S.spm_path,[S.contrasts{c} '_clusters']);
        end
    end

    % For each contrast
    for cldir = 1:length(S.clus_path)

        cdata = dir(fullfile(S.clus_path{cldir},gclusname));

        % For each cluster
        for c = 1%1:length(cdata)
            disp(['*** Cluster = ' num2str(c) ' / ' num2str(length(cdata))])

            [~,cname,~] = fileparts(cdata(c).name);

            % load
            Sc = load(fullfile(S.clus_path{cldir},cdata(c).name));
            
            % update S with data from Sc
            f = fieldnames(Sc.S);
            for i = 1:length(f)
               S.(f{i}) = Sc.S.(f{i})
            end
            
            % Get cluster names
            fnames = fieldnames(S.wf);
            
            % only analyse first one for now
            wf=S.wf.(fnames{1}).wf;
            
            % get wf size to reshape later
            sizewf = size(wf);
            
            if S.ana_singlesub
                wf = wf{S.ana_singlesub}';
            else
                % concatenate over subjects and conditions
                wf = cat(2,wf{:})';
            end

            % remove means over time
            %wf = demean(wf, 1);
            
            L=[];
            start_tol=1e-10;
            tol=start_tol;
            while isempty(L)
                disp(['trying tolerance = ' num2str(tol) ' ...'])
                % get linearly dependent columns and their indices
                % ucol is the order of the magnitude of the columns
                % https://uk.mathworks.com/matlabcentral/answers/108835
                [ld_wf,ld_sub,ucol]=licols(wf,tol);
                %uni_ucol = unique(ucol);
                %for u = uni_ucol'
                %    mu = mean(wf(:,ucol==u));
                %    find(ucol==u,1,'first')
                %end

                % orthogonalise (can also do this on a sample only to get the
                % weight matrix)
                %try
                    [L, ~, ~, W,r,asize] = symmetric_orthogonalise_CAB(ld_wf, 1);
                    disp(['rank = ' num2str(r) ', size = ' num2str(asize)])
                %catch
                if (asize-r)/r > 0.05
                    tol=tol*10;
                else
                    tol=tol+start_tol;
                end
            end
            
            % remove means over time
            %L = demean(L, 1);
            
            % perform (optional) hierarchical cluster analysis
           
            if S.ana_singlesub
                L=L'
            else
                % convert back to subjects/conditions
                L=reshape(L',size(L,2),[],sizewf(1));
                % convert back to cells per subj/cond
                L=squeeze(num2cell(L,[1 2]));
            end
            
            % save 
            S.wf.(fnames{1}).orth_wf = L;
            S.wf.(fnames{1}).ucol = ucol;
            save(fullfile(S.clus_path{cldir},[cname '.mat']),'S');
        end
    end
end