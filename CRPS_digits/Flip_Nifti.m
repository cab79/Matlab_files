clear all
%grplist = [1 2; 29 30]; 
%grplist = [33 34; 31 32]; 
grplist = [39 40; 41 42];
no_cond = 5;
epeaks = [2:3];
no_peaks = length(epeaks);
cdir = pwd;
%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/Unflipped_individuallatency - affected groups separately';
%else
%    filepath = 'S:\Data\CRPS_Digit_Perception_exp1\SPM image files\Unflipped_individuallatency - affected groups separately';
%end
if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/LOR individual';
else
    filepath = 'W:\Data\CRPS_Digit_Perception_exp1\SPM image files\LOR individual';
end
loadsubj

cd(filepath)

for g = 1:size(grplist,2)
    subjects = subjlists(grplist(:,g));
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
                    fname = fullfile(filepath,['spm8_' tmp_nme '_' num2str(epeaks(j)) '_' num2str(i)]);
                    flip_lr([fname '.nii'],[fname 'flipped.nii']);
                end
            end
        end
    end
end

