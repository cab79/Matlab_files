function CORE_condor_run_CRPS_it0

pth = pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

filebase = '';
grp = 'CRPS';
it_in = 0;
CORE_condor_fit_step1(pth,filebase,grp,it_in);
quit