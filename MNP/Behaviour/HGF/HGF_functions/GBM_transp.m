function [pvec, pstruct] = GBM_transp(r, ptrans)
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
idx=r.c_prc.priormusi;

for pn=1:length(nme)
    if r.c_prc.varparam(pn) % if it is a variance parameter
        pvec(idx{pn}) = exp(ptrans(idx{pn}));
    else
        pvec(idx{pn}) = ptrans(idx{pn});
    end
    nme2 = strsplit(nme{pn,1},'log');
    pstruct.(nme2{end}) = pvec(idx{pn});
end

return;
