function [D,S,sim] = MNP_sim_HGF(D,S)

S=MNP_perceptual_models(S);   
Dtemp=D;

if ~isfield(S,'sim_param') || isempty(S.sim_param)
    S=MNP_sim_parameters(S)
else
    S.sim=S.sim_param;
end
    

%sim_param = [mu_0 sa_0 al0 al1 rho ka om eta0 eta1 rb];
sim=[];
for ns = 1:S.numsimrep
    sim{ns}=MNP_HGF(D,S,S.sim);
end

for ns = 1:S.numsimrep
    if length(sim{ns}.y)==length(Dtemp.Sequence.condnum)
        D(ns) = Dtemp;
        D(ns).Output.presstrial = 1:length(sim{ns}.y);
        if length(unique(sim{ns}.y+1))==2 % should be binary, otherwise may be RT
            D(ns).Output.pressbutton = D(ns).Output.Settings.buttonopt(sim{ns}.y+1);
        end
    else
        error('simulated data has the wrong number of trials')
    end
end