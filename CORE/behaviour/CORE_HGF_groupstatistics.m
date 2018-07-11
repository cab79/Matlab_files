function [stats] = CORE_HGF_groupstatistics(T,Vec)
% T is a table containing subject and group information as well as column
% variables for analysis (scalar values per variables and participant).
% V is a structure containing variables (fields) as vectors per participant.

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
    stats.ranksum.(Stfields{stf})= ranksum(gdat{1},gdat{2});
    stats.signtest.(Stfields{stf})= signtest([St(:).(Stfields{stf})]);
    
    pvec(:,stf)=[St(:).(Stfields{stf})]';
end

stats.mvc.param = lda_class(pvec,T.groups,length(T.groups),5);

for vi = 1:length(Vec)

    V=Vec{vi};
    
    % concatenate vector for multivariate classification
    Vfields = fieldnames(V);
    vcat=nan(length(V),length(V(1).(Vfields{1})));
    for v = 1:length(V)
        for vf = 1:length(Vfields)
            vcat(v,:) = V(v).(Vfields{vf})';
        end
    end
    
    stats.mvc.traj(vi) = lda_class(vcat,T.groups,length(T.groups),5);

end