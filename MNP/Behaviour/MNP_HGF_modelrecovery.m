function [S,D,fit,evi] = MNP_HGF_modelrecovery(S,D)

all_models=S.perc_model;

fit=[];
for m=1:length(all_models)
    S.perc_model = all_models(m);
    S=MNP_perceptual_models(S);
    for ns = 1:S.numsimrep
        fit{m}{ns}=MNP_HGF(D(ns),S,[]);
        evi.LME(m,ns)=fit{m}{ns}.optim.LME;
        evi.AIC(m,ns)=fit{m}{ns}.optim.AIC;
        evi.BIC(m,ns)=fit{m}{ns}.optim.BIC;
    end
end
save(fullfile(S.path.prep,['evidence_model' num2str(S.model_sim) '.mat']), 'evi','fit');
