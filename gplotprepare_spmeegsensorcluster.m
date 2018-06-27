function D = gplotprepare_spmeegsensorcluster(S)


%% prepare SPM EEG data
D(1).spm_path=S.spm_path;
D(1).clusdir=S.clusdir;
D(1).facplot=S.facplot;
%S.plotclus;
%S.wavetype; % source or sensor?
%S.wfname; %generic cluster waveform file name
%S.batch; %name of batch .mat file saved from design_batch.m and within same folder
%S.subfactname; 

%%
%dbstop if error

%% plot cluster waveforms
% load waveform data
W=load(fullfile(S.spm_path,S.clusdir,S.wfname));
W=W.S;

% load SPM design and list factors
load(fullfile(S.spm_path,S.batch));
S.fact = {matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fac(:).name};

% identify file indices relating to factor of interest and subject
fact_col = [];
for c = 1:length(S.facplot) % using a loop ensures that fact_col is in correct order for condition labels applied later
    fact_col(c) = find(ismember(S.fact,S.facplot(c)));
end
sub_col = find(ismember(S.fact,S.subfactname));
fact_ind = W.Fm(:,1+fact_col);
sub_ind = W.Fm(:,1+sub_col);
if any(strcmp(S.fact,'Time')) && S.timelev
    time_ind = S.Fm(:,1+find(strcmp(S.fact,'Time')));
    fact_ind = fact_ind(time_ind==S.timelev,:);
    sub_ind = sub_ind(time_ind==S.timelev);
end

if ~isempty(S.plotclus)
    % restrict to selected clusters
    clnames = S.plotclus';
else
    % obtain all cluster names
    clnames = fieldnames(W.wf);
end

% load example cluster image to extract time information
Cnii = load_nii(fullfile(S.spm_path,S.clusdir,[clnames{1} '.nii']));
t_res = Cnii.hdr.hist.srow_z(3);
t_off = Cnii.hdr.hist.srow_z(4);
t_dim = Cnii.hdr.dime.dim(4);
t_last = t_off+t_res*(t_dim-1);
itimes = t_off:t_res:t_last;

% load cluster statistics
load(fullfile(S.spm_path,S.clusdir,'cluster_table.mat'));

range=[];
% for each cluster, 
for cl = 1:length(clnames)
    C=strsplit(clnames{cl},'_');
    cllabel = C{1};
    
    % select cluster data
    wf = W.wf.(clnames{cl});
    % select wf if struct (e.g. source data)
    if isstruct(wf)
        itimes = wf.time*1000;
        wf=wf.wf;
    end
    % identify unique rows of the combination of factors
    [~,Frows,D(cl).WFrows] = unique([fact_ind sub_ind],'rows','stable');
    
    % create new factor indices for unique rows.
    D(cl).Fi = fact_ind(Frows,:);
    D(cl).Si = sub_ind(Frows,:);
    
    % average wf data over non-unique rows
    wff = cell(length(unique(D(cl).WFrows)),1);
    for i = unique(D(cl).WFrows)'
        if size(wf{1},1)>size(wf{1},2) % timepoints in first dimension
            wff{i}=mean(cell2mat(wf(D(cl).WFrows==i)'),2);
        else
            wff{i}=mean(cell2mat(wf(D(cl).WFrows==i)),1)';
        end
    end
    
    D(cl).ptitle = cllabel;
    if size(S.cval,1)>1 && isfield(S,'selectlev')
        D(cl).Fi_ind = find(D(cl).Fi(:,1)==S.selectlev); % indices of Fi (and wf) for each plot
%         cond = P.cval{1}(P.cval{2});
%         condind=nan(1,length(D.cond));
%         for cn = 1:length(cond) 
%             ind = ismember(D.cond,cond{cn});
%             condind(ind) = cn;
%         end
    else
        %[~,D(cl).Fi_ind] = sort(D(cl).Fi);
        D(cl).Fi_ind = 1:length(D(cl).Fi);
    end
    
    % extract cluster statistics
    ct_ind = find(ismember(clustable(:,1),cllabel));
    try
        colind = strcmp(clustable(2,:),'Temporal extent (ms)');
        D(cl).E_val = [min([clustable{ct_ind,colind}]),max([clustable{ct_ind,colind}])];
        D(cl).P_val = [clustable{ct_ind,colind}(1)];
    catch
        D(cl).E_val = [clustable{ct_ind,11:12}];
        D(cl).P_val = [clustable{ct_ind,10}];
    end
    
    % construct condition labels (of last or only factor)
    condlev=S.cval{end,1};
    factind=D(cl).Fi(D(cl).Fi_ind,end);
    D(cl).fi=ismember(factind,S.cval{end,2});
    factind = factind(D(cl).fi);
    D(cl).cond = condlev(factind);

    % set x and y axis data to plot
    D(cl).y=wff(D(cl).Fi_ind);
    D(cl).y=D(cl).y(find(D(cl).fi));
    if ~isempty(S.xlimits)
        xlim = dsearchn(itimes',S.xlimits')';% x values of the selected segment to plot
        D(cl).x = itimes(xlim(1):xlim(2));
        for n = 1:length(D(cl).y)
            D(cl).y{n} = D(cl).y{n}(xlim(1):xlim(2));
        end
    else
        D(cl).x=itimes;
    end

    D(cl).fact_names = S.fact_names;
end
%x=range
end