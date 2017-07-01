clear all;
cd('W:\Data\CRPS_Digit_Perception_exp1\Results\SPM source stats\Regression - indCond\flipped\LORind\RT\subject effect\LOR\130ms');
load cov
clear fnames;
%grplist = [33 34 31 32]; %unflipped
%grplist = [35 36 37 38]; %unflipped
%grplist = [1 2 29 30]; 
grplist = [39 40 41 42]; %unflipped
%grplist = [41]; %unflipped
no_cond = 5;
epeaks = [6];
no_peaks = length(epeaks);

cdir = pwd;

if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/LORflip-ind';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data\CRPS_Digit_Perception_exp1\SPM image files\LORflip-ind';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end

subjects = subjlists(grplist);

%cd(filepath)
%grplist = [33 34 31 32]; %unflipped



Ns=0;
subjnmes = cell(1,1);
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        Ns=Ns+1;
        tmp_nme = subjects{s,1}{s2,1};
        tmp_nme = strrep(tmp_nme, '.left', '_left');
        tmp_nme = strrep(tmp_nme, '.Left', '_left');
        tmp_nme = strrep(tmp_nme, '.right', '_right');
        tmp_nme = strrep(tmp_nme, '.Right', '_right');
        tmp_nme = strrep(tmp_nme, '.flip', '_flip');
        tmp_nme = strrep(tmp_nme, '.Flip', '_flip');
        tmp_nme = strrep(tmp_nme, '.aff', '_aff');
        tmp_nme = strrep(tmp_nme, '.Aff', '_aff');
        tmp_nme = strrep(tmp_nme, '.Unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '.unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '_Left', '_left');
        tmp_nme = strrep(tmp_nme, '_Right', '_right');
        tmp_nme = strrep(tmp_nme, '_Flip', '_flip');
        tmp_nme = strrep(tmp_nme, '_Aff', '_aff');
        tmp_nme = strrep(tmp_nme, '_Unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '.Exp1', '_Exp1');
        
        %tmp_nme = strrep(tmp_nme, 'left', 'left_meansub');
        %tmp_nme = strrep(tmp_nme, 'right', 'right_meansub');
        
        Cs = strsplit(tmp_nme,'_');
         subjnme = Cs{1};
         subjnmes{Ns,1} = subjnme;

         if ~exist(subjnme,'dir')
             mkdir(subjnme);
              copyfile('job_template.mat',fullfile(pwd,subjnme,'job.mat'));
         end
        

        for i = 1:no_cond
            for j = 1:no_peaks
                ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
                if strfind(tmp_nme,'right')
                    fnames{ind,1} = fullfile(filepath,['maspm8_flip_' tmp_nme '_' num2str(epeaks(j)) '_' num2str(i) '.nii']); 
                else
                    fnames{ind,1} = fullfile(filepath,['maspm8_' tmp_nme '_' num2str(epeaks(j)) '_' num2str(i) '.nii']);
                end
            end
        end
    end
end

subjnmes = unique(subjnmes);
for s = 18:length(subjnmes)
    if isempty(cov)
        errordlg('cov cannot be empty');
    end
    %C = [C{:}];   % [EDITED] Most likely this is not useful
    IndexC = strfind(fnames, subjnmes{s});
    Index = find(not(cellfun('isempty', IndexC)));
    cd(fullfile(cdir,subjnmes{s}));
    load job.mat
    matlabbatch{1,1}.spm.stats.factorial_design.dir = {pwd};
    matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = fnames(Index);
    matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.mcov.c = cov(Index);
    save job.mat matlabbatch
    
    spm('defaults','eeg');
    spm_jobman('initcfg');
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(pwd,'spm.mat')};
    spm_jobman('run',matlabbatch);
end
 
for s = 1:length(subjnmes)
    cd(fullfile(cdir,subjnmes{s}));
    load(fullfile(cdir,'job_contrast_template.mat'));
    spm('defaults','eeg');
    spm_jobman('initcfg');
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(pwd,'spm.mat')};
    spm_jobman('run',matlabbatch);
end

fnmes = cell(1,1);
cons = {'con_0001.img'};
nf = 0;
for s = 1:length(subjnmes)
    for c = 1:length(cons)
        nf = nf+1;
        fnmes{nf,1} = fullfile(cdir,subjnmes{s},cons{c});
    end
end
cd(cdir);
load('job.mat');
matlabbatch{1,1}.spm.stats.factorial_design.dir = {pwd};
matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = fnmes;
save job.mat matlabbatch
spm('defaults','eeg');
spm_jobman('initcfg');
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(pwd,'spm.mat')};
spm_jobman('run',matlabbatch);
load(fullfile(cdir,'job_grpcontrast_template.mat'));
matlabbatch{1}.spm.stats.con.spmmat = {fullfile(pwd,'spm.mat')};
spm_jobman('run',matlabbatch);


