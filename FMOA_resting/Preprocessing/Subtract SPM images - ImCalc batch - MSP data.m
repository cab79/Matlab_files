clear all

subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};

Nsub = length(subjects);

for n = 1:Nsub
    
     fnames={
    '_3_1.nii';
    '_3_2.nii';
    '_3_3.nii';
    '_3_4.nii';
    };
    
    subject = subjects(n);
    subject = char(subject);   


    for x = 1:length(fnames)
    fname=char(fnames(x));
    fnames(x)= {['sw_mspm_' subject fname]};
    end
    
    Pfile = char(fnames);
     
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_l_2-1.nii'],'i2-i1')
    Output = spm_imcalc_ui(Pfile,['sw_mspm_' subject '_l_4-3.nii'],'i4-i3')


end
