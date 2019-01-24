function CORE_condor_fit_job
load input.mat;
out = tapas_fitModel_CAB(y, u, prc_model, obs_model, opt_algo, S, 0);
save('output.mat','out');