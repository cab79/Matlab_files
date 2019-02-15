#!/bin/bash

pth=$(pwd)

##add main Matlab functions to parent directory
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_eeg_trial_regression_statistics_condor.m $pth
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_eeg_trial_statistics_condor_compile_outputs.m $pth

## run Matlab step 1
matlab_run CORE_eeg_trial_regression_statistics_condor.m

# copy files to parent dir
#cp $pth/Data/input_files/input*.mat $pth
  
## submit job to condor pool
matlab_submit CORE_condor_EEG_job_run.sub

## run Matlab step 2
matlab_run CORE_eeg_trial_statistics_condor_compile_outputs.m