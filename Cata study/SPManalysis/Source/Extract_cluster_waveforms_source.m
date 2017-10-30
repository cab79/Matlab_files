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
    
    % extract time window of analysis
    Css = strsplit(S.spm_paths{sp,1},'_t')
    TWss = strsplit(Css{2},'_');
    TW = cellfun(@str2double,TWss);
    
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
    
    % select Time level
    if S.timelev && strcmp(S.factlev{S.timelev,1},'Time')
        selt = find(S.Fm(:,1+S.timerow)==S.timelev);
        S.Fm = S.Fm(selt,:);
        S.imglist = S.imglist(selt);
    end
    
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
        
        if isempty(cimages); continue; end

        % For each cluster
        S.wf=struct; % clears this variable for the next iteration
        for c = 1:length(cimages)

            [~,cname,~] = fileparts(cimages(c).name);

            % create empty cell array for waveform data
            S.wf.(cname).wf = cell(length(S.imglist),1); 

            % load cluster image
            %Cnii = load_nii(fullfile(S.clus_path{cldir},cimages(c).name));
            Cnii=spm_vol(fullfile(S.clus_path{cldir},cimages(c).name));
            
            % load VOI
            if S.use_VOI
                voi = load(fullfile(S.clus_path{cldir},['VOI_' cimages(c).name(1:2) '.mat']));
            else
                voi=[];
            end
            
            % for each subjects' source results
            for i = 1:length(usfname)
                disp(['contrast ' num2str(cldir) '/' num2str(length(S.clus_path)) ', cluster ' num2str(c) '/' num2str(length(cimages)) ', image ' num2str(i) '/' num2str(length(usfname))]); % display progress
                sname = fullfile(S.datafile_path,usfname{i});
                
                % extract waveforms
                [Ds,Sw] = SPM_sourcewave_extract(sname, S.spm_paths{sp,2}, Cnii, S.sep_clus,voi,S.use_VOI);
                
                Sw = vertcat(Sw{:}); 
                
                % for each condition for this subject, find corresponding indices of S.imglist
                for t = 1:length(Ds.condlist)
                    try
                        ind = intersect(strmatch(Ds.condlist{t},conname),find(ui2==i));
                    catch
                        ind=[];
                    end
                    if ~isempty(ind)
                        S.wf.(cname).wf{ind} = Sw(:,:,t);
                    end
                end
                
            end
            S.wf.(cname).time = Ds.time;
            
            % perform scaling (to provide equivalent results to SPM source images)
            %if scale_to_VOI
            %    % get time window indices
            %    tw_ind = dsearchn(Ds.time',TW'/1000);
            %    % calculate mean with time window
            %    for i = 1:length(S.wf.(cname).wf)
            %        avwin(i,1) = mean(S.wf.(cname).wf{i}(tw_ind)); 
            %    end
            %    % scaling
            %    sc = voi.Y(selt)./avwin;
            %    % re-scale
            %    for i = 1:length(S.wf.(cname).wf)
            %        S.wf.(cname).wf{i} = S.wf.(cname).wf{i}*sc(i); 
            %    end
            %end
            
            % remove empty cells
            rm_ind = cellfun(@isempty,S.wf.(cname).wf);
            S.wf.(cname).wf(rm_ind)=[];
            
        end
        % save 
        save(fullfile(S.clus_path{cldir},['cluster_data.mat']),'S');

    end
    
end