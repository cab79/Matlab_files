clear all

subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M28';'M35';'M36';'M37';'M38';'M40'};
Nsub = length(subjects);

for n = 1:Nsub
    

     fnames={
    '_1_1.nii';
    '_1_2.nii';
    '_1_3.nii';
    '_1_4.nii';
    };
    
    subject = subjects(n);
    subject = char(subject);   


    for x = 1:length(fnames)
    fname=char(fnames(x));
    fnames(x)= {['sw_spm_base_' subject fname]};
    end
    
    Pfile = char(fnames);
     
    Output = spm_imcalc_ui(Pfile,['sw_substracted_' subject '_l_gamma3_open.nii'],'i1-i3')
    Output = spm_imcalc_ui(Pfile,['sw_substracted_' subject '_l_gamma3_closed.nii'],'i2-i4')

end
