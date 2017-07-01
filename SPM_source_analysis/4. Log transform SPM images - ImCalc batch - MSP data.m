clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_'}; 
Nsub = length(subjects);

for n = 1:Nsub
    
    fpre = 'mspm_';
    
     fnames={
    '_1_1.nii';
    '_1_2.nii';
    '_2_1.nii';
    '_2_2.nii';
    '_3_1.nii';
    '_4_2.nii';
    '_5_1.nii';
    '_6_2.nii';
    };
    
    subject = subjects(n);
    subject = char(subject);   


    for x = 1:length(fnames)
    fname=[fpre subject char(fnames(x))];
  
    Output = spm_imcalc_ui(fname,['log_' fname],'log10(i1)');

    end

end
