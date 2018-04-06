clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
    load([char(subjects(sub)) 'total_data_ICA2.mat']);
    load([subject 'data_matrix_dim2.mat']);
    eval(['data = total_data_ICA2;']);

    if x(3)+x(4) ~= size(data,3); errordlg('trials per condition do not match data size');end
    t1 = 1:x(3);
    t2 = x(3)+1:x(3)+x(4);

    for i=1:size(data,3)
        data(:,:,i)=(detrend(squeeze(data(:,:,i))'))';
    end

    data = blcorrect4(data, 250);

    load epoched_data_template2 %SPM-EEG_template_data
    D.trials = [D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials D.trials];
    D.Fsample = 500;
    D.Nsamples = 2750;

    D.trials = D.trials(1:size(data,3));

    for t = 1:size(data,3)
        if ismember(t,t1)
            D.trials(:,t).events.value = 1;
        elseif ismember(t,t2)
            D.trials(:,t).events.value = 2;
        end
        D.trials(:,t).label = num2str(D.trials(:,t).events.value);
        %D.trials.onset = 0; % onset = time of beginning of trial
        D.trials(:,t).events.time = D.trials(:,t).onset + 4; % time of event
    end

    D.data.fnamedat = ['spm_epoch_' char(subjects(sub)) '.dat'];
    D.fname = ['spm_epoch_' char(subjects(sub)) '.mat'];
    D.data.scale = D.data.scale([1:30 33:64]);
    D.channels = D.channels([1:30 33:64]);
    nchan = length(D.channels);

    nsampl = D.Nsamples;
    datafile = file_array(D.data.fnamedat, [nchan nsampl size(data,3)], 'float32-le');
    D.data.y = datafile;
    D.path = pwd;

    datafile(end,end) = 0;
    fname = D.fname;
    datafile(:, :, :) = data(:,:,:);

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


    clear total_data_ICA2 D


end

