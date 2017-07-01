
%clear all
subject = 'OA17';
load ([subject '_data_epochs.mat']);
load(['b_' subject]);
act = b * total_data;
ib = pinv(b);
%eegplot(act,'winlength',50); 
if ~exist('EEG','var')
    EEG = pop_loadset('eeglab_template.set','\\nasr.man.ac.uk\mhsrss$\snapped\unreplicated\hprg\Chris Brown\Chris Data\Expectancy Study')
end
EEG.setname = subject;
EEG.filename = subject;
%EEG.nbchan = 62;
EEG.pnts = 500;
EEG.srate = 250;
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
