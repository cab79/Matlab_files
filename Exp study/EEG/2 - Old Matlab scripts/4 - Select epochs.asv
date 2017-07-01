clear all
%eegplot(total_data_ICA(:,:,:),'winlength',20,'spacing',100);
subject = 'H5_'; 
filename = 'total_data_ica';
load ([subject filename '.mat']);
EEG = pop_loadset('eeglab_template.set','C:\Documents and Settings\mdmoscab\Desktop\Chris Data\Expectancy Study');
EEG.setname = subject;
EEG.filename = subject;
%EEG.nbchan = 62;
eval(['EEG.data = total_data_ICA;']);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.setname = 'set';
%EEG.nbchan = 62;
ALLEEG(1) = EEG;
CURRENTSET = 1;
EEG.reject.rejmanual = [];
pop_eegplot(EEG,1,1,0);
while isempty(EEG.reject.rejmanual)
    pause(1)
end
reject = find(EEG.reject.rejmanual == 1);