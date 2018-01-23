function h = TactileAWE(h,opt)

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'SimpleSequence';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'1Hz'};
    
    case '1Hz'
        
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
    
    %% EQUIPMENT CONTROL
    % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
    % serial port
    h.Settings.serial = 'COM1';
    %h.Settings.record_EEG='serial';
    %h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
    %h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
    h.Settings.labjack=1; % Use labjack for controlling any equipment?
    h.Settings.labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses GoOne.
    h.Settings.stimcontrol='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.labjack_DACport = 4;
    h.Settings.stimchan = [6]; h.Settings.stimchanforLJ = 1;
    h.Settings.nrchannels = 1; % total number of channels, e.g. on sound card
    % channels on stimulator to use; use differenr rows for different pairs
    % (e.g. for different conditions). If labjack, this can refer to output
    % port(s).
    
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
    h.Settings.trialdur = 4; % if 0, consecutive stimuli will occur with no gap
    % duration of stimulus in seconds
    h.Settings.stimdur = 0; % modified by oddball settings
    % Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.patternmethod = '';
    h.Settings.patternvalue = []; % one per stimdur
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    % Audio sampling rate 
    h.Settings.fs = 96000; % don't change this
    % Binarual beats frequency: creates right ear frequency of f0+df
    %h.Settings.df = 0; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
    % Monaural beats instead? 
    %h.Settings.monaural = 0; 
    % attenuation level in decibels
    %h.Settings.atten = -30; 
    % pitch
    %h.Settings.f0 = 200; % Left ear carrier frequency (pitch)
    %intensity
    h.Settings.inten = 10; % value between 2 and 1000mA for Digitimer DS8R
    % Tactile: number of pulses in a train
    h.Settings.npulses_train = 3; % set to zero to be determined by stimdur
    % Tactile: within-train frequency (Hz)
    h.Settings.p_freq=1; 
    
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
    h.Settings.ntrials = 50;
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
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    
    %% ADAPTIVE
    % adaptive staircase: meanings of the buttonopt
    %h.Settings.adaptive.buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    %h.Settings.adaptive.signalval = [1 2];
    % starting level of adaptive staircase
    %h.Settings.adaptive.startinglevel = max(h.Settings.oddballvalue)-min(h.Settings.oddballvalue); % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100)
    % adapt to omissions of response (not suitable for 2AFC tasks!!)
    %h.Settings.adaptive.omit = 1; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    %h.Settings.adaptive.startomit = 2;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    
    % number of trials of plot
    h.Settings.plottrials=0;
    
    %% THRESHOLDING
    % starting level of adaptive staircase
    h.Settings.threshold.startinglevel = 2; % for intensity
    % adapt to omissions of response (not suitable for 2AFC tasks!!)
    h.Settings.threshold.step = 2;

end

