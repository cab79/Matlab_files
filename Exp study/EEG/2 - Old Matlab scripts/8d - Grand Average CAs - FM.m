subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16'};

subjects = subjects([17:32]);

NSub = length(subjects);



fnames={'_1_avg_ca.mat';
    '_2_avg_ca.mat';
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';
    '_5_avg_ca.mat';
    '_6_avg_ca.mat';
  };

for x=1:length(fnames);
    
grand_avg=zeros(64,2750);

for n = 1:NSub

    fname=char(fnames(x));
    subject=char(subjects(n));
    fname= [subject fname];
    load(fname);
    grand_avg=grand_avg + avg;
    clear avg       
end

grand_avg = grand_avg / NSub; 

eval(['save ' fname(5:6) 'FM_grand_avg.mat'  ' grand_avg'])
clear grand_avg
end



