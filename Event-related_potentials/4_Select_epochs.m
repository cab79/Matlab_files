clear all
subject = 'P35_'; 
filename = 'total_data_ICA';
load ([subject filename '.mat']);
EEG = pop_loadset('S4_3.set','I:\PET study\EEG session\eeglab test');
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
pop_eegplot(ALLEEG,1,1,0);
reject = find(ALLEEG.reject.rejmanual == 1);