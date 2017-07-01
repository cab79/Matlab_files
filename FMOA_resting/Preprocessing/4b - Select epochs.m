subject = 'OA17';
load ([subject '_total_data_ICA_ca.mat']);
if ~exist('EEG','var')
    EEG = pop_loadset('eeglab_template.set','\\nasr.man.ac.uk\mhsrss$\snapped\unreplicated\hprg\Chris Brown\Chris Data\Expectancy Study')
end
EEG.setname = subject;
EEG.filename = subject;
%EEG.nbchan = 62;
EEG.pnts = size(total_data_ICA,2);
EEG.trials = size(total_data_ICA,3);
EEG.data = total_data_ICA;
EEG.srate=250;
EEG.setname = 'set';
%EEG.nbchan = 62;
ALLEEG(1) = EEG;
CURRENTSET = 1;
assignin('base', 'EEG', EEG);
evalin('base','detectArtefactsFinished=0;');
cmd = [ ...
    '[tmprej tmprejE] = eegplot2trial( TMPREJ,' num2str(EEG.pnts) ',' num2str(EEG.trials) ');' ...
    'EEG = pop_rejepoch(EEG, tmprej, 1);' ...
    'detectArtefactsFinished=1;' ...
    ] ;
%Draw the data.
eegplot(EEG.data,...
    'srate',EEG.srate,...
    'events',EEG.event,...
    'command',cmd,...
    'butlabel','Reject',...
    'winlength',20,...
    'spacing',100);
%Wait until the user has finished reviewing.
reviewFinished=0;
while ~reviewFinished
    reviewFinished=evalin('base','detectArtefactsFinished');
    pause(0.01);
end
evalin('base','clear detectArtefactsFinished');
EEG=evalin('base','EEG');
rej = find(tmprej==1);
