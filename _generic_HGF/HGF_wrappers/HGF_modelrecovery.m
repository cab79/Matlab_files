function [S,D,evi] = HGF_modelrecovery(S,D)

all_models=S.perc_model;

for d = 1:length(D)
    fit=[];
    for m=1:length(all_models)
        S.perc_model = all_models(m);
        S=S.perc_model_func(S);
        for ns = 1:S.numsimrep
            S.HGF.selectrep = ns;
            D(d)=HGF_run(D(d),S,2);
            evi.LME(m,ns)=D(d).HGF(ns).fit.optim.LME;
            evi.AIC(m,ns)=D(d).HGF(ns).fit.optim.AIC;
            evi.BIC(m,ns)=D(d).HGF(ns).fit.optim.BIC;
        end
    end
    D(d).fit=fit;
    D(d).evi=evi;
end
% average over subjects
if length(D)>1
    fit = mean(D(:).fit);
    evi = mean(D(:).evi);
end
