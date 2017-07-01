clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
subject = 'S11_'; 
filename = 'total_data_ICA3';
load ([subject filename '.mat']);
EEG = pop_loadset('S4_3.set','C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG from PET sessions\eeglab test');
EEG.setname = subject;
EEG.filename = subject;
EEG.nbchan = 62;
eval(['EEG.data = ' filename ';']);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.setname = 'set';
EEG.nbchan = 62;
ALLEEG(1) = EEG;
CURRENTSET = 1;
EEG.reject.rejmanual = [];
pop_eegplot(EEG,1,1,0);
while isempty(EEG.reject.rejmanual)
    pause(1)
end
reject = find(EEG.reject.rejmanual == 1);