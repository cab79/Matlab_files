#!/bin/bash

pth=$(pwd)

##save file stating which group to run
#echo $1 > grp.txt 

##add main Matlab functions to parent directory
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_condor_EEG_step1.m $pth
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_condor_EEG_step2.m $pth

## run Matlab step 1
matlab_run CORE_condor_fit_step1.m

#nsub=$(<nsub.txt)
#nrand=$(<nrand.txt)

for it in {1..16}
do

  # save iteration number
  echo $it > it.txt

  # copy files to parent dir
  cp $pth/Data/input_files/input*.mat $pth
  
  ## submit job to condor pool
  matlab_submit CORE_condor_fit_job_run.sub
  
  ## run Matlab step 2
  matlab_run CORE_condor_fit_step2.m
done

