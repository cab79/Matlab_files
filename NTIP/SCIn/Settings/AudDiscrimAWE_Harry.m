function h = AudDiscrimAWE_Harry(h,opt)


% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'CreateSequence';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'No_Entrain_VolDiscrim','10Hz_Entrain_VolDiscrim','40Hz_Entrain_VolDiscrim','No_Entrain_PitchDiscrim','10Hz_Entrain_PitchDiscrim','40Hz_Entrain_PitchDiscrim'};
    
    case 'No_Entrain_PitchDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'pitch'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [200 400];
    % set general options
    h = setgeneral(h);
    
    case '10Hz_Entrain_PitchDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'pitch'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [200 400];
    % set general options
    h = setgeneral(h);
    
    case '40Hz_Entrain_PitchDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 40; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'pitch'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [200 400];
    % set general options
    h = setgeneral(h);
    
    case 'No_Entrain_VolDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [0.5 1];
    % set general options
    h = setgeneral(h);
    
    case '10Hz_Entrain_VolDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [0.5 1];
    % set general options
    h = setgeneral(h);
    
    case '40Hz_Entrain_VolDiscrim'
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 40; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [0.5 1];
    % set general options
    h = setgeneral(h);
    

end

function h = setgeneral(h)

%% TRIALS or CONTINUOUS?
h.Settings.design = 'continuous';
% if continuous, how many trials ahead should be in the player schedule?
% (applied to stimulation via soundcard only)
h.Settings.ntrialsahead = 2;  %0 = all trials

%% EXPERIMENTAL CONDIITIONS
% name the settings that define orthogonal condtions at a different row
h.Settings.conds = {};

%% Output options
% save sinwave from all trials as part of stim sequence file
h.Settings.savesinwave = 0;

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
%h.Settings.record_EEG='serial';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
%h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=0; % Use labjack for controlling any equipment?
h.Settings.stimcontrol='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
h.Settings.stimchan = [1 2]; % channels on stimulator to use

%% BLOCKING/RUN OPTIONS
% 'divide' = equally divide trials by nblocks; 
% 'cond' = separate block for each condition
h.Settings.blockopt = 'cond';
% further options for 'divide':
    % number of blocks (containing multiple conditions)
%    h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
    %distribute conditions equally among blocks
%    h.Settings.distblocks = 1;
% options to start sequence at beginning of every run
% 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
% cell array
h.Settings.blockstart = {'buttonpress'}; % audio,labjack,audio
% names of any audiofiles
h.Settings.audiofile = {'instruct.wav'}; % labjack

%% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
% duration of stimulus sequence in seconds
h.Settings.totdur = 240; 
% duration of trial in seconds
h.Settings.trialdur = 0; % if 0, consecutive stimuli will occur with no gap
% duration of stimulus in seconds
h.Settings.stimdur = 0.8; % modified by oddball settings
% Pattern type method: intensity, pitch. Not supported: channel, duration
h.Settings.patternmethod = '';
h.Settings.patternvalue = []; % one per stimdur
% 'rand' or 'reg' spacing?
h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
% sampling rate
h.Settings.fs = 96000; % don't change this
% Monaural beats instead? 
h.Settings.monaural = 0; 
% attenuation level in decibels
h.Settings.atten = -30; 
% pitch
h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
%intensity
h.Settings.inten = 1; % value between 0 and 1

%% Condition-dependent stimulus parameters
% Condition method: intensity, pitch, channel
h.Settings.conditionmethod = {};
h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
% Oddball method: intensity, pitch, channel
h.Settings.oddballtype = 'roving'; % options: 'roving', 'classical'


%% SEQUENCE
% Change probablity (CP): each condition is in rows
h.Settings.oddprob = [
    % standard vs oddball
    0.8 0.2
    ];

% index of row of oddprob that are standards and oddballs
h.Settings.standardind = 1; 
h.Settings.oddind = 2; 
% keep oddball trials apart by at least sep_odd standards
h.Settings.sep_odd = 2;
% for sep_odd, which indices of h.Settings.oddballvalue to consider
% each time? (each list will be separated separately)
h.Settings.sep_odd_ind = {[1 2]};
% for each set, ensure a number of leading standards 
h.Settings.std_lead = 0;
% number of minimum sets to randomised together
h.Settings.n_set = []; % 1 = use min set size
% min number of oddballs within each CP condition
h.Settings.n_odd = [60]; % overrides h.Settings.totdur
% min number of oddballs per randomised set, per CP
h.Settings.n_odd_set = [10]; % overrides h.Settings.totdur
% randomise sets?
h.Settings.rand_set = [1]; 

%% RESPONSE PARAMETERS
% record responses during experiment? 0 or 1
h.Settings.record_response = 1;
% buttonpress options: key: keyboard inputs. Blank for no button press
h.Settings.buttontype='key';
% range of keyboard presses indicating a recordable response
h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 

%% ADAPTIVE: General
% which ones to run? (i.e. indices of h.Settings.adaptive)
h.Settings.adaptive_general.adapttypes = [1];
% alternate or randomise runs over types? Alt must have equal number of
% runs for each adapttype
h.Settings.adaptive_general.seqtype = 'rand'; % 'alt' or 'rand'
h.Settings.adaptive_general.seqrandblocksize = 10; % should divide the number of trials in a set
h.Settings.adaptive_general.selectcond.cp = [1]; % which CP condition to run adaptive on?

%% ADAPTIVE

% adapt on every trial or only just before an oddball?
h.Settings.adaptive(1).oddonly = 1;
% max number of trials after oddball that subject must respond (otherwise counts as omitted response)
h.Settings.adaptive(1).resptrials = 6;


    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 400;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(1);sqrt(sqrt(1))]*2;
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = max(h.Settings.oddballvalue)-min(h.Settings.oddballvalue); % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(1).omit = 1; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(1).startomit = 2;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(1).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(1).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(1).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = max(h.Settings.oddballvalue)-min(h.Settings.oddballvalue); % should be a DIFFERENCE value in mA. 

% number of trials of plot
h.Settings.plottrials=0;
