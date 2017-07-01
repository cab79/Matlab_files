clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
subject = 'S1_'; 
load ([subject 'data_segments2.mat']);
%eegplot(total_data);
load([subject 'b2' ]);
act = b * total_data;
ib = pinv(b);
%eegplot(act,'winlength',50); 
EEG = pop_loadset('S4_3.set','C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG from PET sessions\eeglab test')
EEG.setname = subject;
EEG.filename = subject;
EEG.nbchan = 62;
EEG.pnts = 1501;
EEG.trials = size(total_data,2)/EEG.pnts;
EEG.data = total_data;
EEG.icaweights = b;
EEG.icawinv = ib;
EEG.icaact = act;
EEG.setname = 'set';
EEG.nbchan = 62;
EEG.reject.gcompreject = zeros(30,1);
EEG.stats.compenta = [];
ALLEEG(1) = EEG;
CURRENTSET = 1;
EEG = pop_selectcomps(EEG, [1:30]);
reject = find(EEG.reject.gcompreject==1)'
