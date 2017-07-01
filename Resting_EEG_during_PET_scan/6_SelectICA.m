clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
subject = 'S11_'; 
load ([subject 'data_segments3_acc.mat']);
varname = 'total_data_ICA3';
load([subject 'b3_acc' ]);
eval(['act = b * ' varname ';']);
ib = pinv(b);
%eegplot(act,'winlength',50); 
EEG = pop_loadset('S4_3.set','C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG from PET sessions\eeglab test')
EEG.setname = subject;
EEG.filename = subject;
EEG.nbchan = 62;
EEG.pnts = 1501;
eval(['EEG.trials = size(' varname ',2)/EEG.pnts']);
eval(['EEG.data = ' varname ';']);
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
