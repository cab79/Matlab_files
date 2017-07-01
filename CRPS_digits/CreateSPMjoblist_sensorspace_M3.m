%clear all;
%grplist = [1 3 10 12]; %flipped
clear fnames;
%grplist = [1 2 29 30]; %unflipped left vs right
%grplist = [33 34 31 32]; %unflipped affected vs unaffected
%grplist = [47:50];
grplist = [35:38]; %left right exp1
%grplist = [39 40 41 42]; %unflipped
%grplist = [41]; %unflipped
%usecond = {1, 'dx3', 5};
usecond = {1, 2,3,4, 5};
no_cond=length(usecond);
epeaks = [1];
no_peaks = length(epeaks);
use_flipped = 0;
use_flip_reorient = 1;
suffix = '';

cdir = pwd;


filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\SPM image files\Sensorspace_images';
%filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\SPM image files\eleposnorm\Sensorspace_images';
    run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
%end


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
        
        %tmp_nme = strrep(tmp_nme, 'left', 'left_meansub');
        %tmp_nme = strrep(tmp_nme, 'right', 'right_meansub');
        
        %tmp_nme = ['maspm8_' tmp_nme];
        
        %if strfind(tmp_nme, 'left')
        %    trials = [1:5];
        %elseif strfind(tmp_nme, 'right')
        %    trials = [6:10];
        %end
            
        %for i = 1:no_cond
        %    for j = 1:no_peaks
        %        ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
        %        if use_flipped==1
        %            if strfind(tmp_nme, 'right')
        %                nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*flip_reorient.img'));
        %            else
        %                nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*.img'));
        %            end
        %        else
        %            nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],'strial*.img'));
        %        end
        %        fnames{ind,1} = fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],nme(1).name);
        %    end
        %end
        
        if use_flipped==1
            if strfind(tmp_nme, 'left')
                tmp_nme = ['maspm8_' tmp_nme suffix];
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                tmp_nme = ['maspm8_flip_' tmp_nme suffix];
                trials = [6:10];
            end
        else
            if strfind(tmp_nme, 'left')
                trials = [1:5];
            elseif strfind(tmp_nme, 'right')
                trials = [6:10];
            end
            tmp_nme = ['maspm8_' tmp_nme suffix];
        end
            
        for i = 1:no_cond
            if ischar(usecond{i})
                for j = 1:no_peaks
                    ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
                    fnames{ind,1} = fullfile(filepath,tmp_nme,['smean_' usecond{i} '.img']);
                end
            else
                for j = 1:no_peaks
                    ind = (Ns-1)*no_cond*no_peaks + (j-1)*no_cond + i;
                    if ~isempty(strfind(tmp_nme, 'right')) && use_flip_reorient
                        nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(usecond{i}))],'strial*flip_reorient.img'));
                    else
                        nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(usecond{i}))],'strial*.img'));
                    end
                    
                    fnames{ind,1} = fullfile(filepath,tmp_nme,['type_' num2str(trials(usecond{i}))],nme(1).name);
                end
            end
        end
    end
end
cd(cdir);
load job.mat
matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = fnames;
save job.mat matlabbatch