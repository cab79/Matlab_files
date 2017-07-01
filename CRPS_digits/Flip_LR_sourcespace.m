%clear all;

clear fnames;
%grplist = [33 34 31 32]; %unflipped
grplist = [35 36 37 38]; %unflipped
%grplist = [1 2 29 30]; 
%grplist = [39 40 41 42]; %unflipped
%grplist = [41]; %unflipped
no_cond = 5;
epeaks = [1:4];
no_peaks = length(epeaks);

cdir = pwd;


    filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\SPM image files\eleposnorm';
    %filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\SPM image files\eleposnorm';

    run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');


subjects = subjlists(grplist);

cd(filepath)

Ns=0;
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
        
            
        for i = 1:no_cond
            for j = 1:no_peaks
                ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
                file = dir(fullfile(filepath,['maspm8_' tmp_nme '_' num2str(epeaks(j)) '*' num2str(i) '.nii']));
                fname = fullfile(filepath,file.name); 
                %if strfind(fname, 'right')
                    [pth, nme, ext] = fileparts(fname);
                    flip_lr(fname, fullfile(pth,[nme '_flip.nii']));
                    Output = spm_imcalc_ui({fullfile(pth,[nme '.nii']);fullfile(pth,[nme '_flip.nii'])},fullfile(pth,[nme '_flip_reorient.nii']),'(i1+i2)-i1'); %re-orients image for SPM analysis
                %end
            end
        end
    end
end