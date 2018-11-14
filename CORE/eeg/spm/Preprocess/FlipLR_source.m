clear all

restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA
filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\sourceimages_GS'; 

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of folders
%fpref = 't-200_899_b-200_0_mspm12';
fpref = 'mspm12_fnums';
fmid = '';
%fsuff = '_4_cleaned_tm';
%fsuff = '_2_merged_cleaned';
fsuff = '_4_merged_cleaned'; % fnum

% timewins to flip
tcond = [1 2];

% conditions to flip (right stim only):
%fcond = [5:8 13:16 21:24]; % part4
fcond = [5:8]; % part4 fnum
%fcond = [3 4 7 8 11 12]; % part2

%% RUN
fbase=[fpref '*' fmid  '*' fsuff];
fbase=strrep(fbase,'**','*');

for nt = 1:length(tcond)

    for nf = 1:length(fcond)
        files = dir(fullfile(filepath,[fbase '_' num2str(tcond(nt)) '*' num2str(fcond(nf)) '.nii']));
        for f = 1:length(files)
            fname = fullfile(filepath,files(f).name); 
            [pth, nme, ext] = fileparts(fname);
            sname=fullfile(pth,[nme '_flip.nii']);
            flip_lr(fname, sname);
            spm_imcalc_ui({fname;sname},sname,'(i1+i2)-i1'); %re-orients image for SPM analysis
%             if exist(fname,'file') && exist(sname,'file')
%                 delete(fname)
%             end
        end

    end
end
