function SW_connectivity_results(S)

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
contypes = {'partialCorrelationRegularized','partialCorrelation','correlation'};
ptype = 'FDRp';

%% run

if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

roi1 = [];
roi2 = [];
posneg = [];
tval = [];
pval = [];
fweval = [];
fdrval = [];
connval = [];
firstlevel = [];
grouplevel = [];

% For each analysis (time window)
last_tw=0;
for sp = 1:size(S.spm_paths,1)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    
    % run analysis on combined clusters?
    if strcmp(S.gclusname,'comb_clus.nii')
        timewin = ['Timewin_' num2str(S.spm_paths{sp,2})];
        anapath = fullfile(S.spmstats_path,timewin);
        S.contrasts = {'Combined'};
        % continue if this timewindow was just analysis
        if last_tw==S.spm_paths{sp,2}
            continue
        end
        last_tw = S.spm_paths{sp,2};
    else
        anapath = S.spm_path;
    end
    
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
                    S.(f{i}) = Sc.S.(f{i});
                end
            end
            
            % Get cluster names
            fnames = fieldnames(S.wf);
            
            % only analyse first one for now
            cmats = S.wf.(fnames{1}).correlationMats;
            
            % for each frequency
            nFreq = length(cmats);
            for iFreq = 1:nFreq
                % for each first level contrast
                ncont = length(cmats{iFreq}.groupLevel);
                for ct = 1:ncont
                    % for each connectivity type
                    nconn = length(contypes);
                    for cn = 1:nconn
                        val = cmats{iFreq}.groupLevel(ct).(contypes{cn}).(ptype);
                        tv = cmats{iFreq}.groupLevel(ct).(contypes{cn}).T;
                        pu = cmats{iFreq}.groupLevel(ct).(contypes{cn}).p;
                        fdr = cmats{iFreq}.groupLevel(ct).(contypes{cn}).FDRp;
                        fwe = cmats{iFreq}.groupLevel(ct).(contypes{cn}).FWEp;
                        ngrpcon=size(val,3);
                        ndir=size(val,4);
                        for con = 1:ngrpcon
                            for pn = 1:ndir
                                % find significant p values
                                % identify connected ROIs and direction of the effect
                                % (i.e. for which condition it is stronger)
                                idx=find(val(:,:,con,pn)<0.05);
                                if ~isempty(idx)
                                    
                                    [r1,r2] = ind2sub(size(val(:,:,con,pn)),idx);
                                    roi1 = [roi1;r1];
                                    roi2 = [roi2;r2];
                                    tvi=tv(:,:,con,pn); tval = [tval;tvi(idx)];
                                    pui=pu(:,:,con,pn); pval = [pval;pui(idx)];
                                    fwei=fwe(:,:,con,pn); fweval = [fweval;fwei(idx)];
                                    fdri=fdr(:,:,con,pn); fdrval = [fdrval;fdri(idx)];
                                    
                                    posneg = [posneg;repmat(pn,size(r1))];
                                    connval = [connval;repmat(cn,size(r1))];
                                    firstlevel = [firstlevel;repmat(ct,size(r1))];
                                    grouplevel = [grouplevel;repmat(con,size(r1))];
                                end
                            end
                        end
                    end
                end
            end
            
        end
        % save as excel
        results = table(connval,roi1,roi2,posneg,tval,pval,fweval,fdrval,firstlevel,grouplevel);
        next = datestr(now,30);
        fname = fullfile(S.clus_path{cldir},['Connectivity_run_' timewin '_' next '.xlsx']);
        writetable(results,fname,'Sheet',1)

        % create connectivity table 
        if S.use_aal
            load(fullfile(S.clus_path{cldir},'aal_labels.mat'));
            ROI = genvarname(lab(:,2));
            Thead = table(ROI);
            
            % choose connectiity type
            ana_connval=1;
            ana_ind = find(connval==1);
            
            % identify all analyses
            alllevel = [posneg,firstlevel,grouplevel];
            [urows iA iB] = unique(alllevel(ana_ind),'rows');
            nAna = size(urows,1);
            
            % for each analysis
            for an = 1:nAna
                ind = find(iB==an);
                connmat{an,1} = zeros(size(lab,1));
                for i = ana_ind(ind)'
                    connmat{an}(roi1(i),roi2(i)) = 1;
                end
                connmat{an}(logical(eye(size(connmat{an})))) = nan;
                T = array2table(connmat{an},'VariableNames',genvarname(lab(:,2)));
                T = [Thead T];
                writetable(T,fname,'Sheet',1+an)
            end
            
        end
        
    end
end



