clear all;
%grplist = [1 2 29 30]; %flipped
%clear fnames;
%grplist = [33 34 31 32]; %unflipped Exp2
%grplist = [35 37 36 38]; %unflipped
grplist = [39 40 41 42]; %unflipped Exp1
%grplist = [47:50];
%grplist = [41]; %unflipped
no_cond = 5;
meancond = {[2 3 4]};
epeaks = [1];
no_peaks = length(epeaks);
use_flipped=1;
suffix = '_change';

cdir = pwd;

%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/alltrials/SPM image files/Sensorspace_images';
%    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
%else
    filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\SPM image files\Sensorspace_images';
    run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
%end

subjects = subjlists(grplist);

cd(filepath)

Ns=0;
for s = 3:length(subjects)
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
  
        for mc = 1:length(meancond)
            fnames = [];
            for i = meancond{mc}
                for j = 1:no_peaks
                    ind = (j-1)*no_cond + i;
                    nme = dir(fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],['strial*.img']));
                    for n = 1:length(nme)
                        fnames{length(fnames)+1,1} = fullfile(filepath,tmp_nme,['type_' num2str(trials(i))],nme(n).name);
                    end
                end
            end
            fname_out{mc,1} = fullfile(filepath,tmp_nme,['smean_dx' num2str(length(meancond{mc})) '.img']);
            Output = spm_imcalc_ui(fnames,fname_out{mc,1},'mean(X)',{1,[],[],[]});
        end
    end
end


    