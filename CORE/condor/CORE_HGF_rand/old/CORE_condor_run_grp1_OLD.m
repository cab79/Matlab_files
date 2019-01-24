function CORE_condor_run_grp1

pth = pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

filebase = '';
grp = 1;
it_in = 0;
CORE_condor_fit_step1(pth,filebase,grp,it_in);
quit