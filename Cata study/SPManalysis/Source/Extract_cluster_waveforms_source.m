function Extract_cluster_waveforms_source(S)
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
    
    load(fullfile(S.spm_path,'SPM.mat'));
    S.Fm = SPM.xX.I; % factor matrix
    S.imglist = SPM.xY.P; % Subject-condition image list
    
    % create list of unique spm data file names (usfname)
    % not in same order as in S.imglist
    subname={}; %subject
    conname={}; %trial/condition
    sfname={}; %filename
    for i = 1:length(S.imglist)
        iname=S.imglist{i};
        [~,nme,~] = fileparts(iname);
        C=strsplit(nme,'_');
        subname{i,1} = C{S.subname_index};
        conname{i,1} = C{end};
        sfname{i,1} = [S.prefix subname{i,1} S.spm_paths{sp,3}];
    end
    [usfname ui1 ui2] = unique(sfname,'stable');

    S.clus_path={};
    if isempty(S.contrasts)
        alldir=dir(fullfile(anapath,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(anapath,alldir(c).name);
        end
    else
        for c = 1:length(S.contrasts)
            cona = strrep(S.contrasts{c},' ','');
            cona = strrep(cona,'*','_');
            S.clus_path{c} = fullfile(anapath,[cona '_clusters']);
        end
    end

    % For each contrast
    for cldir = 1:length(S.clus_path)

        cimages = dir(fullfile(S.clus_path{cldir},S.gclusname));

        % For each cluster
        S.wf=struct; % clears this variable for the next iteration
        for c = 1:length(cimages)

            [~,cname,~] = fileparts(cimages(c).name);

            % create empty cell array for waveform data
            S.wf.(cname).wf = cell(length(S.imglist),1); 

            % load cluster image
            %Cnii = load_nii(fullfile(S.clus_path{cldir},cimages(c).name));
            Cnii=spm_vol(fullfile(S.clus_path{cldir},cimages(c).name));
            
            % for each subjects' source results
            for i = 1:length(usfname)
                disp(['contrast ' num2str(cldir) '/' num2str(length(S.clus_path)) ', cluster ' num2str(c) '/' num2str(length(cimages)) ', image ' num2str(i) '/' num2str(length(usfname))]); % display progress
                sname = fullfile(S.datafile_path,usfname{i});
                
                % extract waveforms
                [Ds,Sw] = SPM_sourcewave_extract(sname, S.spm_paths{sp,2}, Cnii, S.sep_clus);
                
                Sw = vertcat(Sw{:}); 
                
                % for each condition for this subject, find corresponding indices of S.imglist
                for t = 1:length(Ds.condlist)
                    try
                        ind = intersect(strmatch(Ds.condlist{t},conname),find(ui2==i));
                    catch
                        ind=[];
                    end
                    if ~isempty(ind)
                        S.wf.(cname).wf{ind(S.timelev)} = Sw(:,:,t);
                    end
                end
            end
            S.wf.(cname).time = Ds.time;
            
            % remove empty cells
            rm_ind = cellfun(@isempty,S.wf.(cname).wf);
            S.wf.(cname).wf(rm_ind)=[];
        end
        % save 
        save(fullfile(S.clus_path{cldir},['cluster_data' num2str(sp) '.mat']),'S');

    end
    
end