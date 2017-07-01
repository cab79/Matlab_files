%clear all;
%grplist = [1 3 10 12]; %flipped
clear fnames;
%grplist = [33 34 31 32]; %unflipped
%grplist = [35 37 36 38]; %unflipped

%grplist = [39 40 41 42]; %unflipped
grplist = [41]; %unflipped
no_cond = 1;
epeaks = [8];
no_peaks = length(epeaks);

cdir = pwd;

if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/LOR individual';
else
    filepath = 'S:\Data\CRPS_Digit_Perception_exp1\SPM image files\LOR individual';
end

loadsubj
subjects = subjlists(grplist);

cd(filepath)
%grplist = [33 34 31 32]; %unflipped

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
        
        tmp_nme = strrep(tmp_nme, 'left', 'mean');
        tmp_nme = strrep(tmp_nme, 'right', 'mean');
        %for i = 1:no_cond
            for j = 1:no_peaks
                ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond +1;
                fnames{ind,1} = fullfile(filepath,['spm8_' tmp_nme '_' num2str(epeaks(j)) '.nii']);
            end
        %end
    end
end
cd(cdir);
load job.mat
matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans = fnames;
save job.mat matlabbatch