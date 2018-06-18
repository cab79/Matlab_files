function [pstruct, pvec] = logrt_softmax_binary_namep(psim,c)

% replace 
pvec = c.priormus;

pstruct = struct;

nme=c.pnames;
nme_mod=c.pnames_mod; % model-specific names
idx=c.priormusi;

for pn=1:length(nme)
    if isfield(psim,nme_mod{pn})
        pstruct.(nme_mod{pn}) = psim.(nme_mod{pn});
        pvec(idx{pn}) = psim.(nme_mod{pn});
    else
        pstruct.(nme_mod{pn}) = pvec(idx{pn});
    end
end

% pstruct.be0 = pvec(1);
% pstruct.be1 = pvec(2);
% pstruct.be2 = pvec(3);
% pstruct.be3 = pvec(4);
% pstruct.be4 = pvec(5);
% pstruct.be5 = pvec(6);
% pstruct.ze  = pvec(7);
% pstruct.be  = pvec(8);

return;
