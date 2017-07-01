clear all
%grplist = [1 2; 29 30]; 
%grplist = [35 36; 37 38]; 
grplist = [39 40; 41 42];
no_cond = 5;
epeaks = [6];
no_peaks = length(epeaks);
cdir = pwd;
%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/Unflipped_individuallatency - affected groups separately';
%else
%    filepath = 'S:\Data\CRPS_Digit_Perception_exp1\SPM image files\Unflipped_individuallatency - affected groups separately';
%end
if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/MSP';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data\CRPS_Digit_Perception_exp1\SPM image files\MSP';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end

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
                    fnames{ind,g} = fullfile(filepath,['maspm8_' tmp_nme '_' num2str(epeaks(j)) '_' num2str(i) '.nii']);
                end
            end
        end
    end
end

for n = 1:size(fnames,1)
  
     files={
    fnames{n,1};
    fnames{n,2};
    };

    Nfile = files{1,1};
    Nfile = strrep(Nfile,'left','mean');
    Nfile = strrep(Nfile,'right','mean');
    Pfile = char(files);
    
    %if ~exist(Nfile,'file');
        Output = spm_imcalc_ui(Pfile,Nfile,'(i1+i2)/2')
    %end
    
    Pfile2 = cell(3,1);
    Pfile2{1,1} = Pfile(1,:);
    Pfile2{2,1} = Pfile(2,:);
    Pfile2{3,1} = Nfile;
    Pfile2 = char(Pfile2);
    
    Nfile1 = files{1,1};
    Nfile1 = strrep(Nfile1,'left','left_meansub');
    Nfile1 = strrep(Nfile1,'right','right_meansub');
    Nfile2 = files{2,1};
    Nfile2 = strrep(Nfile2,'right','right_meansub');
    Nfile2 = strrep(Nfile2,'left','left_meansub');
    
    %if ~exist(Nfile1,'file');
        Output = spm_imcalc_ui(Pfile2,Nfile1,'(i1-i3)')
    %end
    %if ~exist(Nfile2,'file');
        Output = spm_imcalc_ui(Pfile2,Nfile2,'(i2-i3)')
    %end
end