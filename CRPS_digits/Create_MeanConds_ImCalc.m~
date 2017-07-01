clear all
%grplist = [1 2; 29 30]; 
%grplist = [33 34; 31 32]; 
grplist = [39 40; 41 42];
no_cond = 5;
epeaks = [4];
no_peaks = length(epeaks);
cdir = pwd;
%if isunix
%    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/Unflipped_individuallatency - affected groups separately';
%else
%    filepath = 'S:\Data\CRPS_Digit_Perception_exp1\SPM image files\Unflipped_individuallatency - affected groups separately';
%end
if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception_exp1/SPM image files/MSP';
else
    filepath = 'W:\Data\CRPS_Digit_Perception_exp1\SPM image files\MSP';
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
                    fnames{ind,g} = fullfile(filepath,['spm8_' tmp_nme '_' num2str(epeaks(j)) '_' num2str(i) '.nii']);
                end
            end
        end
    end
end

for n = 1:5:size(fnames,1)-4
  fnames_temp = fnames;
     files1={
    fnames_temp{n,1};
    fnames_temp{n,2};
    fnames_temp{n+1,1};
    fnames_temp{n+1,2};
    fnames_temp{n+2,1};
    fnames_temp{n+2,2};
    fnames_temp{n+3,1};
    fnames_temp{n+3,2};
    fnames_temp{n+4,1};
    fnames_temp{n+4,2};
    };

    Nfile = files1{1,1};
    Nfile = strrep(Nfile,'left','mean');
    Nfile = strrep(Nfile,'right','mean');
    Nfile = strrep(Nfile,'_1.nii','.nii');
    Pfile = char(files1);
    
    %if ~exist(Nfile,'file');
        Output = spm_imcalc_ui(Pfile,Nfile,'(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10)/10')
    %end
  
   fnames_temp = fnames;
     files2={
    fnames_temp{n,1};
    fnames_temp{n+1,1};
    fnames_temp{n+2,1};
    fnames_temp{n+3,1};
    fnames_temp{n+4,1};
    };

    Nfile = files2{1,1};
    Nfile = strrep(Nfile,'left','left_mean');
    Nfile = strrep(Nfile,'right','right_mean');
    Nfile = strrep(Nfile,'_1.nii','.nii');
    Pfile = char(files2);
    
    %if ~exist(Nfile,'file');
        Output = spm_imcalc_ui(Pfile,Nfile,'(i1+i2+i3+i4+i5)/5')
    %end
    
    fnames_temp = fnames;
     files3={
    fnames_temp{n,2};
    fnames_temp{n+1,2};
    fnames_temp{n+2,2};
    fnames_temp{n+3,2};
    fnames_temp{n+4,2};
    };

    Nfile = files3{1,1};
    Nfile = strrep(Nfile,'left','left_mean');
    Nfile = strrep(Nfile,'right','right_mean');
    Nfile = strrep(Nfile,'_1.nii','.nii');
    Pfile = char(files3);
    
    %if ~exist(Nfile,'file');
        Output = spm_imcalc_ui(Pfile,Nfile,'(i1+i2+i3+i4+i5)/5')
    %end
end
