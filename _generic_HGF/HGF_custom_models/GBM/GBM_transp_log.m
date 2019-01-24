function [pvec, pstruct] = GBM_transp_log(r, ptrans)
% transforms from native values to log values (only needed for empirical priors script)

pvec    = NaN(1,length(ptrans));
pstruct = struct;

%CAB: names and indices
nme=r.c_prc.pnames;
nme_mod=r.c_prc.pnames_mod; % model-specific names
idx=r.c_prc.priormusi;

for pn=1:length(nme)
    if r.c_prc.varparam(pn) % if it is a variance parameter
        pvec(idx{pn}) = log(ptrans(idx{pn}));
    else
        pvec(idx{pn}) = ptrans(idx{pn});
    end
    pstruct.(nme_mod{pn}) = pvec(idx{pn});
end

return;
