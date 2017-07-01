% CURRENT BASELINE IS LATE ANTICIPATION

clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};

for subject = 1:length(subjects)

d1 = [char(subjects(subject)) 'avg_1_ca'];
d2 = [char(subjects(subject)) 'avg_2_ca'];

data = [];
for i = 1:2
    ii = num2str(i);
    eval(['load(d' ii ')']);
    avg = avg(:,1:2750);
    avg = blcorrect4(avg, 1750); % BASELINE CORRECTED TO PRE-LEP
    data(:,:,i) = avg([1:62],:);
end

load averaged_data_template %SPM-EEG_template_data
D.timeOnset = -4;
D.Fsample = 500;
D.Nsamples = 2750;
D.channels = D.channels(:,[1:30 33:64]);
D.trials = [D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials];
D.trials = D.trials(1:size(data,3));
for notri = 1:size(data,3)
D.trials(1,notri).label = num2str(notri);
end
%sensors = load('H:\Data analysis files\EEG analysis Matlab\LORETA\Electrode coordinates for matlab.txt');
%D.sensors.eeg.pnt = sensors;

direc = 'C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG session';
D.path = direc;

D.data.fnamedat = ['mspm_' char(subjects(subject)) '_LEP.dat'];
D.fname = ['mspm_' char(subjects(subject)) '_LEP.mat'];
nchan = length(D.channels);
nsampl = D.Nsamples;
no_trials = size(data,3);

datafile = file_array(D.data.fnamedat, [nchan nsampl no_trials], 'float32-le');
D.data.y = datafile;
datafile(end,end) = 0;
fname = D.fname;
datafile(:, :, :) = single(data);

D = meeg(D);

S1 = [];
S1.task = 'defaulteegsens';
S1.updatehistory = 0;
S1.D = D;

D = spm_eeg_prep(S1);

save(D)

end




