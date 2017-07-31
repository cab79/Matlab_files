function pstruct = tapas_sdt_namep(pvec)

% CAB run config to get parameter names and indices
C=strsplit(mfilename,'_namep');
eval(['c = ' C{1} '_config;']);

pstruct = struct;

nme=c.pnames;
idx=c.priormusi;

for pn=1:length(nme)
    nme2 = strsplit(nme{pn,1},'log');
    pstruct.(nme2{end}) = pvec(idx{pn});
end

return;
