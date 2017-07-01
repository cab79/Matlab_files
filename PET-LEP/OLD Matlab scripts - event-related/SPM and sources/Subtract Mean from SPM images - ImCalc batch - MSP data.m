clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P14_';'P16_';'P20_';'P22_';'P23_';'P24_';'P25_';'P27_';'P30_';'P31_';'P32_';'P33_';'P35_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S20_'};

Nsub = length(subjects);

for n = 1:Nsub
        
    subject = subjects(n);
    subject = char(subject);   

     fnames={
    '_5_1.nii';
    '_5_2.nii';
    '_6_1.nii';
    '_7_2.nii';
    };

    for x = 1:length(fnames)
    fname=char(fnames(x));
    fnames(x)= {['mspm_' subject '_LEP' fname ]};
    end
    
    Pfile = char(fnames);
     
    Output = spm_imcalc_ui(Pfile,['mspm_' subject '_mean_6-5_1.nii'],'(i1+i3)/2')
    Output = spm_imcalc_ui(Pfile,['mspm_' subject '_mean_7-5_2.nii'],'(i2+i4)/2')
    
     fnames2={
    '_mean_6-5_1.nii';
    '_mean_6-5_2.nii';
    };
    


    for x = 1:length(fnames2)
    fname2=char(fnames2(x));
    fnames2(x)= {['mspm_' subject fname2]};
    end
    
    fnames2 = cat(1,fnames,fnames2);
    
    Pfile2 = char(fnames2);

    Output = spm_imcalc_ui(Pfile2,['mspm_' subject '_meansub_6_1.nii'],'(i3-i5)')
    Output = spm_imcalc_ui(Pfile2,['mspm_' subject '_meansub_7_2.nii'],'(i4-i6)')
    Output = spm_imcalc_ui(Pfile2,['mspm_' subject '_meansub_5_1.nii'],'(i1-i5)')
    Output = spm_imcalc_ui(Pfile2,['mspm_' subject '_meansub_5_2.nii'],'(i2-i6)')
    
    %spm_jobman('initcfg');
    %load smooth_batch
    %matlabbatch{1,1}.spm.spatial.smooth.data = {['mspm_' subject '_meansub_6_1.nii'];['mspm_' subject '_meansub_7_2.nii'];['mspm_' subject '_meansub_5_1.nii'];['mspm_' subject '_meansub_5_2.nii']};
    %matlabbatch{1,1}.spm.spatial.smooth.fwhm = [8 8 8];
    %spm_jobman('run_nogui',matlabbatch);
end
