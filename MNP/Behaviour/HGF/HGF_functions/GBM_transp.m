function [pvec, pstruct] = GBM_transp(r, ptrans)
% transforms from log values back to native values

% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2015 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

pvec    = NaN(1,length(ptrans));
pstruct = struct;

%CAB: names and indices
nme=r.c_prc.pnames;
nme_mod=r.c_prc.pnames_mod; % model-specific names
idx=r.c_prc.priormusi;

for pn=1:length(nme)
    if r.c_prc.varparam(pn) % if it is a variance parameter
        pvec(idx{pn}) = exp(ptrans(idx{pn}));
    else
        pvec(idx{pn}) = ptrans(idx{pn});
    end
    pstruct.(nme_mod{pn}) = pvec(idx{pn});
end

return;
