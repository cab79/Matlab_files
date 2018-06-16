function [pstruct, pvec] = GBM_namep(psim,c)

% replace 
pvec = c.priormus;

pstruct = struct;

nme=c.pnames;
nme_mod=c.pnames_mod; % model-specific names
idx=c.priormusi;

for pn=1:length(nme)
    if isfield(psim,nme_mod{pn})
        pstruct.(nme_mod{pn}) = psim.(nme_mod{pn});
        pvec(idx{pn}) = psim.(nme_mod{pn})(1:length(idx{pn}));
    else
        pstruct.(nme_mod{pn}) = pvec(idx{pn});
    end
end

return;
