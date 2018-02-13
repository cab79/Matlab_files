function h = TENS(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'Ascend'};
    
    case 'Ascend'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'SimpleSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
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
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % one per stimdur
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    % Binarual beats frequency: creates right ear frequency of f0+df
    %h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    % Monaural beats instead? 
    %h.Settings.monaural = 0; 
    % attenuation level in decibels
    %h.Settings.atten = -30; 
    % pitch
    %h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
    %intensity
    h.Settings.inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.maxinten = 200;
    %intensity
    h.Settings.inten_diff = []; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.t_freq=[]; 
    
   
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [];
    h.Settings.oddballtype = ''; % options: 'roving', 'classical'

    
    %% SEQUENCE
    h.Settings.ntrials = 500;
    % Change probablity (CP): each condition is in rows
    %h.Settings.oddprob = [
        % standard (left) vs oddball (right)=
    %    ];
    % index of row of oddprob that are standards and oddballs
    %h.Settings.standardind = 1; % does not apply to 'roving oddball' design
    % keep oddball trials apart by at least sep_odd standards
    %h.Settings.sep_odd = [2 1]; % for each CP condition
    % for each set, ensure a number of leading standards 
    %h.Settings.std_lead = [0 0]; % for each CP condition
    % number of minimum sets to randomise together
    %h.Settings.n_set = [1 1]; % 1 = use min set size, for each CP condition
    % min number of oddballs within each CP condition
    %h.Settings.n_odd = [20, 20]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    %h.Settings.n_odd_set = [5, 5]; % overrides h.Settings.totdur
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'UpArrow','DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    h.Settings.threshold.startinglevel = 0; % for intensity)
    h.Settings.threshold.step = 2;
    h.Settings.threshold.signalval = [2 1]; % 2 = carrying on increasing; 1 = decrease
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
%h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=0; % Use labjack for controlling any equipment?
h.Settings.stimcontrol='serial'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
h.Settings.spt_port = 'COM3';
h.Settings.spt1_port = 'COM5'; % only needed for 'serial' EEG triggers
%h.Settings.labjack_DACport = 4;
%h.Settings.stimchan = [6]; h.Settings.stimchanforLJ = 1;
h.Settings.nrchannels = 1; % total number of channels, e.g. on sound card
% channels on stimulator to use; use differenr rows for different pairs
% (e.g. for different conditions). If labjack, this can refer to output
% port(s).

% Tactile: number of pulses in a train
h.Settings.npulses_train = 1; % set to zero to be determined by stimdur
% Tactile: within-train frequency (Hz)
h.Settings.p_freq=[];
% Use timer to control frequency of labjack outputs? Otherwise uses software timing.
%h.Settings.labjack_timer=1; % must use if there are a large number of pulses (e.g. >4)

% Audio sampling rate 
h.Settings.fs = 96000; % don't change this