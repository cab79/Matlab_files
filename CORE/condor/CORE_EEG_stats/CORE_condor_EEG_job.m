function CORE_condor_EEG_job
load input.mat;
PG=S.pred_group;
[out,S] = bayesreg_crossval(X,Y,S,PG);
save('output.mat','out','S','X','PG','stats');