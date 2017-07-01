clear all;
%grplist = [1 3 10 12]; %flipped
%clear fnames;
%grplist = [33 34 31 32]; %unflipped Exp2
%grplist = [35 37 36 38]; %unflipped
grplist = [39 40 41 42]; %unflipped Exp1
%grplist = [47:50];
%grplist = [41]; %unflipped
no_cond = 5;
meancond = {[1 5]; [2 3 4]};
epeaks = [1];
no_peaks = length(epeaks);

cdir = pwd;


filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials';
run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m')


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
       
        
        tmp_nme = ['aspm8_' tmp_nme];
        trials = [1:5];

            
        for mc = 1:length(meancond)
            fnamesall = [];
            fc=0;
            for i = meancond{mc}
                fc = fc+1;
                for j = 1:no_peaks
                    ind = (j-1)*no_cond + i;
                    nme = dir(fullfile(filepath,[tmp_nme '*.nii']));
                    % find number of unique files for each condition number
                    indC=(strfind({nme.name},[num2str(trials(i)) '.nii']));
                    ind = find(not(cellfun('isempty', indC)));
                    % add into fnames
                    for f = 1:length(ind)
                        fnamesall{fc,f} = fullfile(filepath,nme(ind(f)).name);
                    end
                end
            end
            for f=1:size(fnamesall,2)
                fnames = fnamesall(:,f);
                fname_out{mc,1} = strrep(fnames{1},'.nii',['mean_c' num2str(mc) '.nii']);
                Output = spm_imcalc_ui(fnames,fname_out{mc,1},'mean(X)',{1,[],[],[]});
            end
        end
        
    end
end


    