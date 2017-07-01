cd('M:\Matlab\ExpStudy');
clear all

rmv_out=100; % values that are a factor of rmv_out * median are removed from trajectories

ranked =0;
max_stim = 360;
rname = 'ExpStudy_HGF_results_BU_v2.mat';
tname = 'ExpStudy_HGF_results_BU_v2.xlsx';
load(rname);
results = orderfields(results);

subs = fieldnames(results);
nsub =length(subs);

num_stim = 2;
var = {'al','alpf','dau','da','ud','wt','psi','epsi','mu','sa','AIC','BIC','LME','irr'};
fol = {'p_prc','p_prc','traj','traj','traj','traj','traj','traj','traj','traj','optim','optim','optim','irr'};
numvar = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %2 level model
%numvar = [1, 1, 3, 1, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1]; %3 level model
num_stimtypes = [num_stim*ones(1,8) ones(1,6)];
num_cond = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, 6 = extract each condition
incl_abs = [0 0 1 1 1 0 0 1 0 0 0 0 0 0];
nme_list = {};
sub_sort=[16,24,25,26,27,28,29,30,31,17,18,19,20,21,22,23,1,9,10,11,12,13,14,15,2,3,4,5,6,7,8,32,41,42,43,44,45,46,47,33,34,35,36,37,38,39,40];

mu_traj = nan(nsub,max(numvar),num_stim*max_stim);

for v=1:length(var)
    for nv = 1:numvar(v)
        for ns = 1:num_stimtypes(v)
            if numvar(v)>1
                if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns) '_' num2str(nv)];
                    eval([nme ' = nan(nsub,1);']);
                else
                    nme = [var{v} num2str(nv)];
                    eval([nme ' = nan(nsub,1);']);
                end
            else
                 if num_stimtypes(v)>1
                    nme = [var{v} num2str(ns)];
                    eval([nme ' = nan(nsub,1);']);
                 else
                    nme = var{v};
                    eval([nme ' = nan(nsub,1);']);
                 end
            end
            nme_list{length(nme_list)+1} = [nme '_av']; 
            if incl_abs(v)==1
                nme_list{length(nme_list)+1} = [nme '_abs']; 
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


for s = 1:nsub;
    sub = subs{s};
    cond = results.(sub).cond;
    stimlen = length(results.(sub).traj.dau);
    inds = [];
    for i = 1:num_stim
        inds(i,:) = i:num_stim:(stimlen-num_stim+i);
    end
    
    for nv = 1:max(numvar)
        a = results.(sub).traj.mu(:,nv);
        mu_traj(s,nv,1:length(a)) = a;
    end
    mu_traj=mu_traj(sub_sort,:,:);
    
    for v=1:length(var)
        for nv = 1:numvar(v)
            for ns = 1:num_stimtypes(v)
                if numvar(v)>1
                    if num_stimtypes(v)>1
                        nme = [var{v} num2str(ns) '_' num2str(nv)];
                        if strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results.(sub).(fol{v}).' [var{v} num2str(ns)] '(nv);']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results.(sub).(fol{v}).' var{v} '(inds(ns,:),nv);']);
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
                        if strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results.(sub).(fol{v}).' var{v} '(nv);']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results.(sub).(fol{v}).' var{v} '(:,nv);']);
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
                        if strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results.(sub).(fol{v}).' [var{v} num2str(ns)] ';']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results.(sub).(fol{v}).' var{v} '(inds(ns,:));']);
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
                            eval(['a = results.(sub).' var{v} '(inds(ns,:));']);
                            %a(a>abs(rmv_out*median(a)))=[];
                            %a(a<-abs(rmv_out*median(a)))=[];
                            eval([nme '(s,1) = length(a);']);
                        end 
                     else
                        nme = var{v};
                        if strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
                            eval([nme '(s,1) = results.(sub).(fol{v}).' var{v} ';']);
                        elseif strcmp(fol{v},'traj');
                            eval(['a = results.(sub).(fol{v}).' var{v} ';']);
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
                            eval(['a = results.(sub).' var{v} ';']);
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
T.Subs = subs;
for nv = 1:length(nme_list)
    eval(['T.' nme_list{nv} '=' nme_list{nv};]);
end
T=T(sub_sort,:);
writetable(T,tname);
