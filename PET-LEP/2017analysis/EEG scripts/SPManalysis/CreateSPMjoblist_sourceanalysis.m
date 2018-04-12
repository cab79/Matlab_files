clear fnames;
filepath = 'C:\Data\PET-LEP\SPM_source_images\LOR_ind_20170222T224400';
anapath = 'C:\Data\PET-LEP\SPM_source_images';
run('C:\Data\Matlab\Matlab_files\PET-LEP\EEG scripts\loadsubj.m');
grplist = [1,2,3,... % healthy
            4,5,6];   % patient
no_cond = 2;
epeaks = [1 2 10 11 12 13];
no_peaks = length(epeaks);
use_flip_reorient = 0;

cdir = pwd;


subjects = subjlists(grplist);

cd(filepath)

Ns=0;
for s = 1:length(subjects)
    for s2 = 1:size(subjects{s,1},1) 
        Ns=Ns+1;
        tmp_nme = subjects{s,1}{s2,1};
        for i = 1:no_cond
            for j = 1:no_peaks
                ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
                %if ~isempty(strfind(tmp_nme,'right')) && use_flip_reorient
                %    file = dir(fullfile(filepath,['maspm8_' tmp_nme '_' num2str(epeaks(j)) '*' num2str(i) '_flip_reorient.nii']));
                %    fnames{ind,1} = fullfile(filepath,file.name); 
                %else
                    file = dir(fullfile(filepath,['spm12_' tmp_nme '_' num2str(epeaks(j)) '_*' num2str(i) '.nii']));
                    fnames{ind,1} = fullfile(filepath,file.name);
                %end
            end
        end
    end
end
cd(cdir);
load job.mat
matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = fnames;
save job.mat matlabbatch