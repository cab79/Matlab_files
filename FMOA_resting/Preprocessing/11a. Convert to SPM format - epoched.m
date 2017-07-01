clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA3','OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16'};


for subject = 1:length(subjects)

d1 = [char(subjects(subject)) '_total_data_ICA_ca.mat'];
d2 = [char(subjects(subject)) '_block_info.mat'];
load(d1); load(d2); 

clear total_data_ICA

for i=1:size(total_data_ICA,3)
    total_data_ICA(:,:,i)=detrend(squeeze(total_data_ICA(:,:,i))')';
end

load epoched_data_template2 %SPM-EEG_template_data
D.trials = [D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials];

ii = find(events_mat(:,2)== 1 & events_mat(:,3)==1 | events_mat(:,2)==2 & events_mat(:,3)==1);

data = zeros(size(total_data_ICA,1),size(total_data_ICA,2),length(ii));
for i = 1:length(ii)
data(:,:,i) = total_data_ICA(:,:,ii(i));
end

clear total_data_ICA

D.trials = D.trials(1:size(data,3));
D.Nsamples = 1000;
D.timeOnset = 0;

for t = 1:length(ii)
D.trials(:,t).events.value = events_mat(ii(t),2);
D.trials(:,t).label = num2str(D.trials(:,t).events.value);
%D.trials.onset = 0; % onset = time of beginning of trial
D.trials(:,t).events.time = D.trials(:,t).onset + 0; % time of event
end

D.data.fnamedat = ['spm_base_' char(subjects(subject)) '.dat'];
D.fname = ['spm_base_' char(subjects(subject)) '.mat'];
D.data.scale = D.data.scale([1:2 4:30 33:64]);
nchan = length(D.channels);
nsampl = D.Nsamples;
datafile = file_array(D.data.fnamedat, [nchan nsampl length(ii)], 'float32-le');
D.data.y = datafile;
D.path = pwd;

datafile(end,end) = 0;
fname = D.fname;
datafile(:, :, :) = data([1:2 4:30 33:64],:,:);

D = meeg(D);
S1 = [];
S1.task = 'defaulttype';
S1.D = D;
S1.updatehistory = 0;
D = spm_eeg_prep(S1);

        %S1 = [];
        %S1.task = 'defaulteegsens';
        %S1.updatehistory = 0;
        %S1.D = D;

        D = spm_eeg_prep(S1);

save(D)


clear ii events_mat data datafile



end

