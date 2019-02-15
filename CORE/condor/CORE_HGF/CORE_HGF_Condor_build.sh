#!/bin/bash

pth=$(pwd)

##unzip
unzip dependencies.zip
#unzip Data.zip

#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_HGF_rand/CORE_condor_fit_job.m $pth
#cp $pth/dependencies/Matlab_files/CORE/condor/CORE_HGF_rand/CORE_condor_monitor_build.m $pth

## build exe
matlab_build CORE_condor_fit_job.m

## monitor build
matlab_run CORE_condor_monitor_build.m

## clean up
rm -i dependencies.zip
rm -i Data.zip