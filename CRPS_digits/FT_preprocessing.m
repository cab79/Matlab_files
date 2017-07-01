%% This is the EEGLAB

% choose a file

% 'C:\Users\jay\Desktop\Work\Fieldtrip Example data\jimher_toolkit_demo_dataset_.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\Pilot2_pairedpulse_LeftM1_45.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\Pilot2_singlepulse_LeftM1_45.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Giovanni test 03-06-2014\singlepulse_52.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Giovanni test 03-06-2014\lici100.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_80.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_100.vhdr';
% 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_140.vhdr';
addpath('C:\Users\jay\Desktop\Work\Programs');
addpath('C:\Program Files\MATLAB\R2011a\toolbox\eeglab13_3_2b');
eeglab

EEG = pop_loadbv('C:\Users\jay\Desktop\Work\Fieldtrip Example data','jimher_toolkit_demo_dataset_.vhdr');
EEG.setname='single100';
EEG = eeg_checkset( EEG );

%choose channel file
%EEG=pop_chanedit(EEG, 'load',{'C:\\Users\\jay\\Desktop\\Work\\Channels\\jay_good_channels_30.ced' 'filetype' 'autodetect'});
EEG=pop_chanedit(EEG, 'load',{'C:\\Users\\jay\\Desktop\\Work\\Channels\\jay_good_channels.ced' 'filetype' 'autodetect'},'delete',32);
EEG = eeg_checkset( EEG );
makeEvent(EEG,5000);
EEG = pop_importevent( EEG, 'append','no','event','C:\\Users\\jay\\Desktop\\Work\\Programs\\event.txt','fields',{'latency' 'type'},'skipline',1,'timeunit',1);
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  'p-pulse'  }, [-1  1], 'newname', 'single epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-1000     -1]);
EEG = eeg_checkset( EEG );


%%    First make sure that you use EEGLAB to import the events and epoch the 
%     data. When that EEG file is in the workspace this should work
% % %
%real shit
addpath('C:\Users\jay\Desktop\Work\fieldtrip-20140804');
addpath('C:\Users\jay\Desktop\Work\fieldtrip-20140804\fileio');
%% first we load the cfg with the data and define the trial
%
%  ** make sure to change the dataset when you use different data sets
cfg = [];
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\Pilot2_pairedpulse_LeftM1_45.vhdr';
cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\Pilot2_singlepulse_LeftM1_45.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Giovanni test 03-06-2014\singlepulse_52.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\Fieldtrip Example data\jimher_toolkit_demo_dataset_.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Giovanni test 03-06-2014\lici100.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_80.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_100.vhdr';
%cfg.dataset = 'C:\Users\jay\Desktop\Work\EEG Tests\Second Practice test 02-07-2014\ramp_leftM!_140.vhdr';

cfg.hdr = ft_read_header(cfg.dataset);
cfg.continuous              = 'no';
cfg.trialdef.prestim        = 1;         % prior to event onset
cfg.trialdef.poststim       = 1;        % after event onset
cfg.trialdef.eventtype      = 'p-pulse'; % see above
cfg.trialdef.eventvalue     = 1000 ;

%     Change this when not doing field trip example
cfg.data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
data = cfg.data;
cfg.trl = ft_makeEvent(cfg);
%cfg.trialfun = 'ft_makeEvent';    %might need to use this if the above line doesnt work

cfg = ft_definetrial(cfg);
trl = cfg.trl;

%% REMOVE BAD CHANNELS
% 
%  ** not that these channels have to be changed. Probably easiest to see
%  in eeglab
selchan = ft_channelselection({'all' '-T7' '-TP10' '-FC1' '-EMG1' '-EMG2'}, cfg.data.label);
data = ft_selectdata(data, 'channel', selchan);

%% Filter 60 Hz noise
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.9 -0.01];
cfg.preproc.bsfilter = 'yes';
cfg.preproc.bsfreq = [59 61];

data = ft_preprocessing(cfg, data);

%% visually reject bad trials
%



cfg =[];
cfg.method = 'trial';
cfg.alim = 1e-5;

data = ft_rejectvisual(cfg, data);

%%
%update trl
 trl = update_trl( data, trl );
% Apply original structure to segmented data, gaps will be filled with nans
cfg     = [];
cfg.trl = trl;
data = ft_redefinetrial(cfg, data); 


%% set rejection markers

cfg = [];
cfg.trl = trl;
cfg.continuous = 'no';
cfg.method = 'marker';
cfg.prestim = 0.000;
[cutoff time] = ft_getCutoff(data, 1)  % 1 if single pulse, 2 if double pulse
cfg.poststim = cutoff;
cfg.Fs = 5000;
cfg.trialdef.eventtype  = 'p-pulse';
cfg.trialdef.eventvalue = 1000;
cfg.trialfun = 'ft_markpulse';
[cfg_artifact, artifact] = ft_artifact_tms(cfg, data);


%% reject artifact
 cfg_artifact.artfctdef.reject = 'partial';
 cfg_artifact.artfctdef.minaccepttim = 0.01;
 data = ft_rejectartifact(cfg_artifact, data); 
 
%% Browse the data 
 
cfg = [];
 cfg.preproc.demean = 'yes';
 cfg.preproc.baselinewindow = [-1 -0.01];
 cfg.layout = 'easycapM23.lay';
 ft_databrowser(cfg, data);
 
 
%% Display the segmented data including the artifacts that are gone
 %
 
 cfg = [];
cfg.artfctdef = cfg_artifact.artfctdef; % Store previously obtained artifact definition
cfg.continuous = 'yes'; % Setting this to yes forces ft_databrowser to represent our segmented data as one continuous signal
ft_databrowser(cfg, data);

% %% visually reject bad trials
% %
% 
% %  *** make sure to reject 2 trials here since they are split in 2
% 
% cfg =[];
% cfg.method = 'trial';
% cfg.alim = 1e-5;
% 
% data = ft_rejectvisual(cfg, data);
% %%
% %update trl
%  trl2 = update_trl( data, data.cfg.previous.trl );
% % Apply original structure to segmented data, gaps will be filled with nans
% cfg     = [];
% cfg.trl = trl2;
% data = ft_redefinetrial(cfg, data); 
%% visually reject LAST TWO TRIAL
%

data = remove_nan_trial(data);

%% Perform ICA on segmented data
%

cfg = [];
cfg.demean = 'yes'; 
cfg.method = 'fastica';        % FieldTrip supports multiple ways to perform ICA, 'fastica' is one of them.
cfg.fastica.approach = 'symm'; % All components will be estimated simultaneously.
cfg.fastica.g = 'gauss'; 
 
comp_tms = ft_componentanalysis(cfg, data);

%% Show averages of the time analysis and topographic images

cfg = [];
cfg.vartrllength  = 2; % This is necessary as our trials are in fact segments of our original trials. This option tells the function to reconstruct the original trials based on the sample-information stored in the data
comp_tms_avg = ft_timelockanalysis(cfg, comp_tms);

figure;
cfg = [];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, comp_tms_avg);

figure;
cfg           = [];
cfg.component = [1:size(comp_tms.label,1)];
cfg.comment   = 'no';
cfg.layout    = 'easycapM23'; % If you use a function that requires plotting of topographical information you need to supply the function with the location of your channels
ft_topoplotIC(cfg, comp_tms);

figure;
cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, comp_tms);

%% Use unmixing matrix to get original data back to remove components

cfg          = [];
cfg.demean   = 'no'; % This has to be explicitly stated as the default is to demean.
cfg.unmixing = comp_tms.unmixing; % Supply the matrix necessay to 'unmix' the channel-series data into components
cfg.topolabel = comp_tms.topolabel; % Supply the original channel label information
 
comp_tms         = ft_componentanalysis(cfg, data);  % MAKE SURE THIS IS SUPPOSED TO BE DATA_VISUAL AND NOT COMP_TMS

%% Reject components

cfg            = [];
cfg.component  = [1 6 18 10 17 20];  
cfg.demean     = 'no';

data_tms_clean_segmented = ft_rejectcomponent(cfg, comp_tms);


%% Restructure the data

% Apply original structure to segmented data, gaps will be filled with nans
cfg     = [];
cfg.trl = trl;
data_tms_clean = ft_redefinetrial(cfg, data_tms_clean_segmented); % Restructure cleaned data


%*** check this out. it might go here or in the next one


%data_tms_clean = remove_nan_trial(data_tms_clean);

%% Reject the last trial
 cfg.preproc.demean = 'yes';
 cfg.preproc.baselinewindow = [-1 -0.01];
 cfg.layout = 'easycapM23.lay';
 ft_databrowser(cfg, data_tms_clean);   
 
 % Reject the last trial

cfg =[];
cfg.method = 'trial';
cfg.alim = 1e-5;

data_tms_clean = ft_rejectvisual(cfg, data_tms_clean);

%% Interpolate the data


% Replacing muscle artifact with nans
    muscle_window = [0.0 cutoff]; % The window we would like to replace with nans
    muscle_window_idx = [nearest(data_tms_clean.time{1},muscle_window(1)) nearest(data_tms_clean.time{1},muscle_window(2))]; % Find the indices in the time vector corresponding to our window of interest
    for i=1:size(data_tms_clean.trial,2) % Loop through all trials
      data_tms_clean.trial{i}(:,muscle_window_idx(1):muscle_window_idx(2))=nan; % Replace the segment of data corresponding to our window of interest with nans
    end;

    % Interpolate nans using cubic interpolation
    cfg = [];
    cfg.method = 'cubic'; % Here you can specify any method that is supported by interp1: 'nearest','linear','spline','pchip','cubic','v5cubic'
    cfg.prewindow = 0.010; % Window prior to segment to use data points for interpolation
    cfg.postwindow = 0.010; % Window after segment to use data points for interpolation
    data_tms_clean = ft_interpolatenan(cfg, data_tms_clean); % Clean data




    %% Filter the data


cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.9 -0.01];
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 80;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.3;
data_filt = ft_preprocessing(cfg, data_tms_clean);

cfg = [];
 cfg.preproc.demean = 'yes';
cfg.preproc.detrend = 'yes';
 cfg.preproc.baselinewindow = [-1 -0.01];
 cfg.layout = 'easycapM23.lay';
 ft_databrowser(cfg, data_filt);
 


    %% Reject bad Epochs

cfg =[];
cfg.method = 'trial';
cfg.alim = 1e-5;
data_filt = ft_rejectvisual(cfg, data_filt);

trl = update_trl( data_filt, trl );
%% ICA again

cfg = [];
cfg.demean = 'yes'; 
cfg.method = 'fastica';        % FieldTrip supports multiple ways to perform ICA, 'fastica' is one of them.
cfg.fastica.approach = 'symm'; % All components will be estimated simultaneously.
cfg.fastica.g = 'gauss'; 
 
comp_tms = ft_componentanalysis(cfg, data_filt);

%% show the components
    cfg = [];
    cfg.vartrllength  = 2; % This is necessary as our trials are in fact segments of our original trials. This option tells the function to reconstruct the original trials based on the sample-information stored in the data
    comp_tms_avg = ft_timelockanalysis(cfg, comp_tms);

    figure;
    cfg = [];
    cfg.viewmode = 'butterfly';
    ft_databrowser(cfg, comp_tms_avg);

    
    cfg           = [];
    cfg.component = [1:size(comp_tms.label,1)];
    cfg.comment   = 'no';
    cfg.layout    = 'easycapM10'; % If you use a function that requires plotting of topographical information you need to supply the function with the location of your channels
    ft_topoplotIC(cfg, comp_tms);
    

cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, comp_tms);
    

%% unmixing 
cfg          = [];
cfg.demean   = 'no'; % This has to be explicitly stated as the default is to demean.
cfg.unmixing = comp_tms.unmixing; % Supply the matrix necessay to 'unmix' the channel-series data into components
cfg.topolabel = comp_tms.topolabel; % Supply the original channel label information
 
comp_tms_new   = ft_componentanalysis(cfg, data_filt);

%% Reject components

cfg            = [];
cfg.component  = [2 15 29 30 31 32 54 24 61 ];
cfg.demean     = 'no' ;
 
data_tms_clean_segmented = ft_rejectcomponent(cfg, comp_tms_new);

%% Restructure the data

%update trl
 trl = update_trl( data_filt, trl );

% Apply original structure to segmented data, gaps will be filled with nans
cfg     = [];
cfg.trl = trl;
data_tms_clean = ft_redefinetrial(cfg, data_tms_clean_segmented); % Restructure cleaned data


 cfg.preproc.demean = 'no'
 cfg.layout = 'easycapM23.lay';
 ft_databrowser(cfg, data_tms_clean);
 


%% Get rid of any nans from the redefine

data_tms_clean = remove_nan_trial( data_tms_clean);


%% Convert back to EEGLAB

fieldtrip2eeglab

%% Remove any nan trials

EEG = remove_nan_trial(EEG);

