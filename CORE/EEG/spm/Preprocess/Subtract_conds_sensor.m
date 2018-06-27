clear all

%% SPECIFY DATA
filepath = 'C:\Data\CORE\SPMdata\sensorimages'; 

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of folders
fpref = 't-200_299_b-200_0_mspm12';
fmid = '';
fsuff = '_4_cleaned_tm';

% conditions to subtract: column 2 minus column 1
subcond = [
    1,3;
    2,4;
    5,7;
    6,8;
    9,11; 
    10,12; 
    13,15;
    14,16;
    17,19;
    18,20;
    21,23;
    22,24];

use_flip=1;

%% RUN
files = dir(fullfile(filepath,[fpref '*' fmid  '*' fsuff]));

for f = 34%1:length(files)
    fname = files(f).name;
 
    % subtract the images
    for nf = 1:size(subcond,1)
        
        if use_flip==1 && exist(fullfile(filepath,files(f).name,['scondition_' num2str(subcond(nf,1)) '_flip.nii']),'file')
            suff= '_flip';
        else
            suff='';
        end
        
        i1 = fullfile(filepath,files(f).name,['scondition_' num2str(subcond(nf,1)) suff '.nii']);
        i2 = fullfile(filepath,files(f).name,['scondition_' num2str(subcond(nf,2)) suff '.nii']);
        outname = fullfile(filepath,files(f).name,['scondition_' num2str(subcond(nf,1)) '-' num2str(subcond(nf,2)) suff '.nii']);
        
        spm_imcalc_ui({i1;i2},outname,'i1-i2');
       
    end
end
