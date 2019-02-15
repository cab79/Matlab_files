function [stats] = CORE_HGF_groupstatistics(T,Vec,S)
% T is a table containing subject and group information as well as column
% variables for analysis (scalar values per variables and participant).
% Vec: cell array. Use trajectories from Vec(1) with names defined from
% Vec(2). If empty, uses traj defined from S.condmean and S.cond

% group indices
grps = unique(T.groups,'stable');
for g = 1:length(grps)
    grpind{g} = find(strcmp(T.groups,grps{g}));
end

% group stats: only works for two groups
St = table2struct(T);
Stfields = fieldnames(St);
ignorefields = {'subjects','groups','sessions','blocks','conds'};
Stfields = Stfields(~ismember(Stfields,ignorefields));
pvec=nan(length(St),length(Stfields));
for stf = 1:length(Stfields)
    if any(isnan([St(:).(Stfields{stf})])) || isempty([St(:).(Stfields{stf})])
        continue
    end
    for g = 1:length(grps)
        for d = 1:length(grpind{g})
            gdat{g}(d) = St(grpind{g}(d)).(Stfields{stf});
        end
    end
    [stats.ranksum.p.(Stfields{stf}),~,stst] = ranksum(gdat{1},gdat{2},'method','approximate');
    stats.ranksum.z.(Stfields{stf}) = stst.zval;
    stats.signtest.(Stfields{stf})= signtest([St(:).(Stfields{stf})]);
    
    pvec(:,stf)=[St(:).(Stfields{stf})]';
end
if isfield(S,'lda') && any(strcmp(S.lda,'para'))
    % remove nans
    pvec = reshape(pvec(~isnan(pvec)),size(pvec,1),{});
    % remove constants
    pvec = pvec(:,~var(pvec,[],1)==0);
    out = lda_class(pvec,T.groups,length(T.groups),5,0);
    try
        stats.mvc.param.error = out.ldaCVErr;

        if S.nperm
            perm_results = lda_perm(S,1:size(pvec),'params',pvec,T,out,stats);
            stats.mvc.param.p_value = sum(abs(perm_results) < abs(stats.mvc.param.error)) / length(perm_results);
        else
            stats.mvc.param.p_value = [];
        end

    catch
        stats.mvc.param.error = [];
        stats.mvc.param.p_value = [];
    end
end

%for vi = 1:length(Vec)

if ~isempty(Vec)
    % use trajectories from Vec(1) with names defined from Vec(2)
    if length(Vec)>1
        Vfields_all = Vec{2};
    else
        Vfields_all = fieldnames(Vec{1});
    end
    Vec = Vec(1);
else
    % use traj defined by S.condmean and S.cond
    conds = fieldnames(S.cond);
    Vfields_all = {};
    for cm = 1:length(S.condmean)
        Vfields_all{cm} = S.condmean{cm};
        vcell={};
        for cd = 1:length(conds)
            vcell = [vcell,{[S.condmean{cm} '_' conds{cd}]}];
        end
        for stf=1:length(vcell)
            for sub = 1:length(St)
                Vec{1}(sub).(Vfields_all{cm})(stf) = St(sub).(vcell{stf});
            end
        end
    end
end
%for cm = 1:length(Vfields_all)

if isfield(S,'lda') && any(strcmp(S.lda,'traj_conds'))
    fn = fieldnames(S.cond);
    Vfields=Vfields_all;
    for vf = 1:length(Vfields)
        dat=[];
        for cond = 1:length(fn)
            dat(:,cond) = T.([Vfields{vf} '_' fn{cond}]);
        end
        out = lda_class(dat,T.groups,length(T.groups),5,0);
        stats.mvc.traj.error.(Vfields{vf}) = out.ldaCVErr;
    end
end

if isfield(S,'lda') && any(strcmp(S.lda,'traj_trials'))
    try
        Vfields=Vfields_all;
        V=Vec{1};
        for vf = 1:length(Vfields)
            dat=[];
            for sub = 1:length(V)
                dat(sub,:) = V(sub).(Vfields{vf})';
            end
            out = lda_class(dat,T.groups,length(T.groups),5,0);
            try
                stats.mvc.traj.error.(Vfields{vf}) = out.ldaCVErr;

                if S.nperm
                    perm_results = lda_perm(S,V,Vfields{vf},dat,T,out, stats);
                    stats.mvc.traj.p_value.(Vfields{vf}) = sum(abs(perm_results) < abs(stats.mvc.traj.error.(Vfields{vf}))) / length(perm_results);
                else
                    stats.mvc.traj.p_value.(Vfields{vf}) = [];
                end

            catch
                stats.mvc.traj.error.(Vfields{vf}) = [];
                stats.mvc.traj.p_value.(Vfields{vf}) = [];
            end
            
            % univariate stats
            for i = 1:size(dat,2)
                stats.ranksum.([Vfields{vf} '_trials'])(i)= ranksum(dat(strcmp(T.groups,grps{1}),i),dat(strcmp(T.groups,grps{2}),i));
            end
            stats.ranksum.([Vfields{vf} '_sigtrials']) = find(stats.ranksum.([Vfields{vf} '_trials'])<0.05);
            
        end
    catch
        return
    end
end
end

function perm_results = lda_perm(S,V,testname,dat,T,in, stats)

perm_results=[];
for np = 1:S.nperm
    ix = randperm(length(V));
    out = lda_class(dat(ix,:),T.groups,length(T.groups),5, 1);
    perm_results(np) = out.ldaCVErr;
    disp(['LDA: ' testname ', perm ' num2str(np) ', error = ' num2str(perm_results(np))])
end

end