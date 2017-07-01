clear all
close all
cd C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\LW\ANOVA

%% settings
configuration = LW_ANOVA('default',[]);
configuration.parameters.factor_data=struct('label',{'group','side','digit'},'within',{0,1,1});
configuration.parameters.permutation=0;
configuration.parameters.enable_parallel = 0;

%% load data and run ANOVA
load datasets % first need to create a 'dataset' structure using the LW GUI
load table_data % first need to save a cell array for the factor matrix
configuration.parameters.table_data=table_data;
update_pointers = [];
[out_configuration,out_datasets] = LW_ANOVA('process',configuration,datasets,update_pointers);

%% save results
dt = datestr(now,30);
save(['LW_ANOVA_' dt],'out_configuration','out_datasets');