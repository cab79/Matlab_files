
if hgf
    fnames = dir(fullfile(aname,['*_bopars' opt_name '.mat']));
    clear rs % struct 
    clear rc % cell
    for f=1:length(fnames)
        load(fullfile(aname,fnames(f).name));
        rs(f)=bopars;
        rc{f}=bopars;
    end
    tapas_fit_plotCorr_group(rs);
    rc=tapas_bayesian_parameter_average_CAB(1,rc);
    tapas_fit_plotCorr(rc);
end