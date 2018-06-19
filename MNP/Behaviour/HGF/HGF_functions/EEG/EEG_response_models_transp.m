function [pvec, pstruct] = EEG_response_models_transp(r, ptrans)
% transforms from log values back to native values

pvec    = NaN(1,length(ptrans));
pstruct = struct;

%CAB: names and indices
nme=r.c_obs.pnames;
nme_mod=r.c_obs.pnames_mod; % model-specific names
idx=r.c_obs.priormusi;

for pn=1:length(nme)
    if r.c_obs.varparam(pn) % if it is a variance parameter
        pvec(idx{pn}) = exp(ptrans(idx{pn}));
    else
        pvec(idx{pn}) = ptrans(idx{pn});
    end
    pstruct.(nme_mod{pn}) = pvec(idx{pn});
end

return;
