subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16'};

subjects = subjects([17:29 31:32]);

NSub = length(subjects);



fnames={'_1_avg_ca.mat';
    '_2_avg_ca.mat';
    '_3_avg_ca.mat';
    '_4_avg_ca.mat';
    '_5_avg_ca.mat';
    '_6_avg_ca.mat';
  };

for x=1:length(fnames);
    
grand_avg=zeros(64,750);

for n = 1:NSub

    fname=char(fnames(x));
    subject=char(subjects(n));
    fname= [subject fname];
    load(fname);
    avg = avg(:,1:750);
    grand_avg=grand_avg + avg;
    clear avg       
end

grand_avg = grand_avg / NSub; 

eval(['save ' fname(5:6) 'FM_grand_avg.mat'  ' grand_avg'])
clear grand_avg
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%late SPN
load('6_HE_grand_avg.mat');
topoplot(mean(grand_avg([1:2 4:30 33:64], 432:500),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

topoplot(mean(CER([1:2 4:30 33:64], 432:500),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'); caxis([-5 5]), colorbar

