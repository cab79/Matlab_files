clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13'};
  
Nsub = length(subjects);

for n = 1:Nsub
    
     fnames={
    '_4_1.nii';
    '_5_2.nii';
    '_6_3.nii';
    '_7_4.nii';
    '_8_5.nii';
    '_9_6.nii';
    '_2_1.nii';
    '_2_2.nii';
    '_2_3.nii';
    '_2_4.nii';
    '_2_5.nii';
    '_2_6.nii';
    '_3_1.nii';
    '_3_2.nii';
    '_3_3.nii';
    '_3_4.nii';
    '_3_5.nii';
    '_3_6.nii';
    };
    
    subject = subjects(n);
    subject = char(subject);   


    for x = 1:length(fnames)
    fname=char(fnames(x));
    fnames(x)= {['sw_mspm_' subject fname]};
    end
    
    Pfile = char(fnames);
     
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_p2_unc.nii'],'(i1+i2+i3)/3')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_p2_cer.nii'],'(i4+i5+i6)/3')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_e_unc.nii'],'(i7+i8+i9)/3')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_e_cer.nii'],'(i10+i11+i12)/3')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_l_unc.nii'],'(i13+i14+i15)/3')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_l_cer.nii'],'(i16+i17+i18)/3')


end
