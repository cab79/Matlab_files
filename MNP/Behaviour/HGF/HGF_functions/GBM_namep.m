function pstruct = GBM_namep(pvec,varargin)

% CAB run config to get parameter names and indices
if nargin>1
    c = varargin{1};
    eval(['c = ' c ';']);
else
    C=strsplit(mfilename,'_namep');
    eval(['c = ' C{1} '_config;']);
end

pstruct = struct;

nme=c.pnames;
idx=c.priormusi;

for pn=1:length(nme)
    nme2 = strsplit(nme{pn,1},'log');
    pstruct.(horzcat(nme2{:})) = pvec(idx{pn});
end

return;
