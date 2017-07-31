clear all
dname='C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmtrials_nan'; num_lev=1;

cd(dname);
dbstop if error
factors = {'Side', 'DC'}; % oly need to label those that 
levels = {{'L' 'R'}, {'1' '3'}}; % oly need to label those that 
levels_aff = {'aff' 'unaff'}; analyse_aff=1;
factor_matrix = [ % must have the same number of rows as there are conditions
          1 1
          1 2
          2 1
          2 2
          1 1
          1 2
          2 1
          2 2
          1 1
          1 2
          2 1
          2 2
          ];
% condition index: same numbers = same condition
[~,~,ci] = unique(factor_matrix,'rows');

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\CORE\Participant_data.xlsx';
% names of headers in the above xls file:
    subhead = 'Subject';
    grphead = 'Group';
    inchead = 'Include';
    CRPSsidehead = 'CRPSlr';
% which codes to analyse in 'Include' columns in participant data file?
include_codes = [1];
[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),grphead));
sub_col = find(strcmp(pdata(1,:),subhead));
inc_col = find(strcmp(pdata(1,:),inchead));
side_col = find(strcmp(pdata(1,:),CRPSsidehead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));

rmv_out=100; % values that are a factor of rmv_out * median are removed from trajectories

ranked =0;
if ranked==1
    tname = 'results_ranked.xlsx';
else
    tname = 'results.xlsx';
end

files_ana = inc_idx;
fi=0;
for f = files_ana' 
    fi=fi+1;
    C = strsplit(pdata{f+1,sub_col},'CORE');
    file = fullfile(dname,[C{2} '_bopars_aff.mat']);
    subs{1,fi} = pdata{f+1,sub_col};
    load(file);
    results(fi)=bopars;
end
results = orderfields(results);
nsub =length(subs);

num_stim = 1; % number of event markers per trial; e.g. anticipation cue AND laser stimulus = 2 stim
ncond = max(ci); % two hands * two DC

var = {'al0','al1','rb','om','dau','da','ud','wt','psi','epsi','mu','sa','AIC','BIC','LME','irr','be0','be1','be2','be3','be4','be5','ze','be'};
fol = {'p_prc','p_prc','p_prc','p_prc','traj','traj','traj','traj','traj','traj','traj','traj','optim','optim','optim','irr','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs'};
num_stimtypes = [num_stim*ones(1,10) ones(1,6) ones(1,8)];
% 3lev
if num_lev==3
    numvar = [1, 1, 0, 3, 1, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; %3 level model
    num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
    incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
elseif num_lev==2
% 2lev, e.g. SDT
    numvar = [1, 1, 1, 0, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0];
    num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
    incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
elseif num_lev==1
% 1lev, e.g. KF
    numvar = [1, 1, 0, 1, 1, 0, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; %3 level model
    num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
    incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

nme_list = {};
sub_sort=1:nsub; % sort order of subjects
grp = pdata(1+files_ana',grp_col);

% remove variables
var(numvar==0)=[];
fol(numvar==0)=[];
num_stimtypes(numvar==0)=[];
num_cond(numvar==0)=[];
incl_abs(numvar==0)=[];
numvar(numvar==0)=[];

mu_traj = nan(nsub,max(numvar),length(results(1).u));

% create variables and list of variable names
V=struct;
for v=1:length(var)
    for nv = 1:numvar(v)
        for ns = 1:num_stimtypes(v)
            if numvar(v)>1
                if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns) '_' num2str(nv)];
                else
                    nme = [var{v} num2str(nv)];
                end
            else
                 if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns)];
                 else
                    nme = var{v};
                 end
            end
            eval([nme ' = nan(nsub,1);']);
            %V.(nme) = nan(nsub,1);
            nme_list{length(nme_list)+1} = [nme '_av']; 
            %V.(nme).name = [nme '_av'];
            if incl_abs(v)==1
                nme_list{length(nme_list)+1} = [nme '_abs']; 
                %V.(nme).name = [nme '_abs'];
            end

            if num_cond(v)>1
                for nc = 1:num_cond(v)
                    nme_list{length(nme_list)+1} = [nme '_' num2str(nc)];
                end
            end
            if incl_abs(v)==1
                if num_cond(v)>1
                    for nc = 1:num_cond(v)
                        nme_list{length(nme_list)+1} = [nme '_' num2str(nc) '_abs']; 
                    end
                end
            end
        end
    end
end

exp = nan(nsub,3);
for s = 1:nsub;
    sub = subs{s};
    cond=[];
    if isfield(results,'conds')
        for c = 1:ncond
            cii = find(ci==c);
            for ciii = 1:length(cii)
                ciiii = results(s).conds==cii(ciii);
                cond(ciiii) = c;
            end
        end
    end
    stimlen = length(results(s).traj.mu);
    
    % indices of each stimulus type (e.g. anticipation cue and pain
    % stimulus)
    inds = [];
    for i = 1:num_stim
        inds(i,:) = i:num_stim:(stimlen-num_stim+i);
    end
    
    for nv = 1:numvar(strcmp(var,'mu'))
        a = results(s).traj.mu(:,nv);
        mu_traj(s,nv,1:length(a)) = a;
    end
    mu_traj=mu_traj(sub_sort,:,:);
    
    for v=1:length(var) % for each variable
        for nv = 1:numvar(v) % for level/type of that variable
            for ns = 1:num_stimtypes(v) % and each stimulus type
                if numvar(v)>1
                    if num_stimtypes(v)>1
                        nme = [var{v} num2str(ns) '_' num2str(nv)];
                        if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results(s).(fol{v}).' [var{v} num2str(ns)] '(nv);']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results(s).(fol{v}).' var{v} '(inds(ns,:),nv);']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            if num_cond(v)>1
                                for nc = 1:num_cond(v)
                                    eval([nme '(s,nc) = nanmean(a(cond==nc));']);
                                end
                            else
                                eval([nme '(s,1) = nanmean(a,1);']);
                            end
                        end
                    else
                        nme = [var{v} num2str(nv)];
                        if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results(s).(fol{v}).' var{v} '(nv);']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results(s).(fol{v}).' var{v} '(:,nv);']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            if num_cond(v)>1
                                for nc = 1:num_cond(v)
                                    eval([nme '(s,nc) = nanmean(a(cond==nc));']);
                                end
                            else
                                eval([nme '(s,1) = nanmean(a,1);']);
                            end
                        end
                    end
                else
                     if num_stimtypes(v)>1
                        nme = [var{v} num2str(ns)];
                        if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results(s).(fol{v}).' [var{v} num2str(ns)] ';']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results(s).(fol{v}).' var{v} '(inds(ns,:));']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            if num_cond(v)>1
                                for nc = 1:num_cond(v)
                                    eval([nme '(s,nc) = nanmean(a(cond==nc));']);
                                end
                            else
                                eval([nme '(s,1) = nanmean(a,1);']);
                            end
                        elseif strcmp(fol{v},'irr');
                            eval(['a = results(s).' var{v} '(inds(ns,:));']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            eval([nme '(s,1) = length(a);']);
                        end 
                     else
                        nme = var{v};
                        if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results(s).(fol{v}).' var{v} ';']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results(s).(fol{v}).' var{v} ';']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            if num_cond(v)>1
                                for nc = 1:num_cond(v)
                                    eval([nme '(s,nc) = nanmean(a(cond==nc));']);
                                end
                            else
                                eval([nme '(s,1) = nanmean(a,1);']);
                            end
                        elseif strcmp(fol{v},'irr');
                            eval(['a = results(s).' var{v} ';']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            eval([nme '(s,1) = length(a);']);
                        end 
                     end
                end
            end
        end
    end

end

for v=1:length(var)
    for nv = 1:numvar(v)
        for ns = 1:num_stimtypes(v)
            if numvar(v)>1
                if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns) '_' num2str(nv)];
                else
                    nme = [var{v} num2str(nv)];
                end
            else
                 if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns)];
                 else
                    nme = var{v};
                 end
            end
            
            % average across conditions
            eval([nme '_av = squeeze(nanmean(' nme ',2));']);
            
            % average within-condition
            if num_cond(v)>1
                for nc = 1:num_cond(v)
                    eval([nme '_' num2str(nc) '= squeeze(nanmean(' nme '(:,nc),2));']);
                end
            end

            % create absolute values
            if incl_abs(v)==1
                eval([nme '_abs = abs(' nme '_av);']);
                if num_cond(v)>1
                    for nc = 1:num_cond(v)
                        eval([nme '_' num2str(nc) '_abs = abs(' nme '_' num2str(nc) ');']);
                    end
                end
            end
            
            %create tied ranks of data
            if ranked == 1;
                % on averages
                eval(['vsize = size(' nme '_av);']);
                eval(['isort = floor(tiedrank(' nme '_av(:)));']);
                eval([nme '_av= reshape(isort,vsize(1),vsize(2));']);
                % on abs averaged data
                if incl_abs(v)==1
                    eval(['vsize = size(' nme '_abs);']);
                    eval(['isort = floor(tiedrank(' nme '_abs(:)));']);
                    eval([nme '_abs= reshape(isort,vsize(1),vsize(2));']);
                end
                
                % on conditions
                if num_cond(v)>1
                    if incl_abs(v)==1
                        eval([nme '_cond_abs = abs(' nme ');']);
                        eval(['vsize = size(' nme '_cond_abs);']);
                        eval(['isort = floor(tiedrank(' nme '_cond_abs(:)));']);
                        eval([nme '_cond_abs= reshape(isort,vsize(1),vsize(2));']);    
                    end
                    eval(['vsize = size(' nme ');']);
                    eval(['isort = floor(tiedrank(' nme '(:)));']);
                    eval([nme '= reshape(isort,vsize(1),vsize(2));']);   
                    for nc = 1:num_cond(v) 
                        if incl_abs(v)==1
                            eval([nme '_' num2str(nc) '_abs= squeeze(nanmean(' nme '_cond_abs(:,nc),2));']);     
                        end
                        eval([nme '_' num2str(nc) '= squeeze(nanmean(' nme '(:,nc),2));']);
                    end
                end
            end
        end
    end
end


T = table;
T.Subs = subs';
T.Grp = grp;
[~,~,T.GrpNum] = unique(grp);
%T.ExpLow = exp(:,1);
%T.ExpMed = exp(:,2);
%T.ExpHigh = exp(:,3);
for nv = 1:length(nme_list)
    eval(['T.' nme_list{nv} '=' nme_list{nv};]);
end
T=T(sub_sort,:);
writetable(T,tname);

% group comparison stats
Ts = table;
for nv = 1:length(nme_list)
    eval(['x=' nme_list{nv} '(T.GrpNum==1);']);
    eval(['y=' nme_list{nv} '(T.GrpNum==2);']);
    if all(~isnan(x))
        eval(['Ts.' nme_list{nv} '= ranksum(x,y);']);
    end
end
writetable(Ts,['Grp_nonpara_' tname]);

% t-test vs zero
Tt = table;
for nv = 1:length(nme_list)
    eval(['x=' nme_list{nv} ';']);
    if all(~isnan(x)) && length(unique(x))>1
        eval(['Tt.' nme_list{nv} '= signtest(x);']);
    end
end
writetable(Tt,['Signtest_' tname]);

