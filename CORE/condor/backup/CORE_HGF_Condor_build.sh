#!/bin/bash

pth=$(pwd)

##unzip
unzip dependencies.zip
unzip Data.zip

cp $pth/dependencies/Matlab_files/CORE/behaviour/CORE_condor_fit_job.m $pth
cp $pth/dependencies/Matlab_files/CORE/behaviour/CORE_condor_monitor_build.m $pth

## build exe
matlab_build CORE_condor_fit_job.m

## monitor build
matlab_run CORE_condor_monitor_build.m

## clean up
rm dependencies.zip
rm Data.zip