clear all
subject = 'H11';
load ([subject '_data_epochs_500Hz.mat']);
load(['b' subject]);
act = b * total_data;
ib = pinv(b);
eegplot(act,'winlength',50); 
EEG = pop_loadset('eeglab_template.set','C:\Documents and Settings\mdmoscab\Desktop\Chris Data\Expectancy Study')
EEG.setname = subject;
EEG.filename = subject;
%EEG.nbchan = 62;
EEG.pnts = 2750;
EEG.trials = size(total_data,2)/EEG.pnts;
EEG.data = total_data;
EEG.icaweights = b;
EEG.icawinv = ib;
EEG.icaact = act;
EEG.setname = 'set';
%EEG.nbchan = 62;
EEG.reject.gcompreject = zeros(30,1);
EEG.stats.compenta = [];
ALLEEG(1) = EEG;
CURRENTSET = 1;
EEG = pop_selectcomps(EEG, [1:30]);
reject = find(EEG.reject.gcompreject==1)'
