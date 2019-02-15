#!/bin/bash

pth=$(pwd)

##unzip
unzip dependencies.zip
#unzip Data.zip

#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_condor_EEG_job.m $pth
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_EEG_stats/CORE_condor_monitor_build.m $pth

## build exe
matlab_build CORE_condor_EEG_job.m

## monitor build
matlab_run CORE_condor_monitor_build.m

## clean up
rm -i dependencies.zip
#rm -i Data.zip