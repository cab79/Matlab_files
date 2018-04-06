clear all
subject = 'P35_'; 
load ([subject 'data_segments2.mat']);
load([subject 'b2' ]);
act = b * total_data_ICA;
ib = pinv(b);
eegplot(act,'winlength',50); 

EEG = pop_loadset('S4_3.set','I:\PET study\EEG session\eeglab test')
EEG.setname = subject;
EEG.filename = subject;
EEG.nbchan = 62;
EEG.pnts = 2750;
EEG.trials = size(total_data_ICA,2)/EEG.pnts;
EEG.data = total_data_ICA;
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

%eegplot(act,'winlength',100); 

reject = find(EEG.reject.gcompreject==1)'
