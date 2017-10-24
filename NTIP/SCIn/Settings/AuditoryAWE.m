function h = AuditoryAWE(h,opt)

%% NOTES ON SETTINGS:
% 1/h.Settings.df*0.25 AND (1/df+1/pitchdiff)*0.25 needs to be a multiple of 1/f0 to prevent clicking caused by phase re-setting of the f0 sound
% Settings that work:
% df=10Hz: f0=200, pitchdiff=200 
% df=25Hz: f0=200, pitchdiff=200 
% df=40Hz: f0=160, pitchdiff=160
%strategy: decide on df first, then pitchdiff (ideally a multiple of df),
%then calculate:
% 1./((1/h.Settings.df+1/pitchdiff)*0.25 ./ [1:10])
% and choose a frequency for f0.
% seems to work if f0 and pitchdiff are both multiples of df, and f0 and
% pitchdiff are also multiples of each other?

% REVISED: f0 and f0+pitchdiff must be a multiple of: 
%    1/((1/h.Settings.df)*0.25) 
%    (for these particular settings - depends on stimdur)

%if ~exist('opt','var')
%    opt = {'10Hz'};
%end

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'CreateSequence';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'10Hz_oddball','1Hz_oddball','10Hz','1Hz'};
    
    case '10Hz'
        h.Settings.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
        h.Settings.patternvalue = [500 500]; % one per stimdur
        h = AuditoryAWE(h,'commonsettings');
        
    case '1Hz'
        h.Settings.df = 1; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
        h.Settings.patternvalue = [500 500]; % one per stimdur
        h = AuditoryAWE(h,'commonsettings');
    
    case '10Hz_oddball'
        h.Settings.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
        h.Settings.patternvalue = [500 540]; % one per stimdur
        h = AuditoryAWE(h,'commonsettings');
        
    case '1Hz_oddball'
        h.Settings.df = 1; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
        h.Settings.patternvalue = [500 540]; % one per stimdur
        h = AuditoryAWE(h,'commonsettings');
        
    case 'commonsettings'    
        %% TRIALS or CONTINUOUS?
        h.Settings.design = 'continuous';
        % if continuous, how many trials ahead should be in the player schedule?
        % (applies to stimulation via soundcard only)
        h.Settings.ntrialsahead = 2;  %0 = all trials

        %% Output options
        % save sinwave from all trials as part of stim sequence file
        h.Settings.savesinwave = 0;

        %% EQUIPMENT CONTROL
        % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
        h.Settings.record_EEG='serial';
        h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
        h.Settings.EEGMarkPattern = 1; % mark EEG for every change in stimulus pattern (0 = start of trial only)
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
        h.Settings.blockstart = {'buttonpress'}; 
        % names of any audiofiles
        h.Settings.audiofile = {''}; 

        %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
        % duration of stimulus sequence in seconds
        h.Settings.totdur = 300; 
        % duration of trial in seconds
        h.Settings.trialdur = 0; % if 0, consecutive stimuli will occur with no gap
        % duration of stimulus in seconds
        h.Settings.stimdur = 0.4; % modified by oddball settings
        % Pattern type method: intensity, pitch. Not supported: channel, duration
        h.Settings.patternmethod = 'pitch';
        % 'rand' or 'reg' spacing?
        h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
        % sampling rate
        h.Settings.fs = 96000; % don't change this
        % Monaural beats instead? 
        h.Settings.monaural = 0; 
        % attenuation level in decibels
        h.Settings.atten = 0; 
        % pitch
        h.Settings.f0 = 500; % Left ear carrier frequency (pitch) - overridden by h.Settings.patternvalue
        %intensity
        h.Settings.inten = 1; % value between 0 and 1

        %% Condition-dependent stimulus parameters
        % Condition method: intensity, pitch, channel
        h.Settings.conditionmethod = {};
        h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
        % Oddball method: intensity, pitch, channel
        h.Settings.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
        % Odball value: DURATION DEVIANTS
        % i.e. all possible duration options and their probability
        % In general, each row is a different stim type, with columns providing
        % within-stimulus temporal patterns of pitch, intensity or channel.
        % Here the options are chosen as:
            % left column = 1st inten/pitch
            % right column = 2nd inten/pitch
        % and the temporal pattern is defined by fc (from either fpitch or finten)
        sd = h.Settings.stimdur;
        h.Settings.oddballvalue = [
            % standard
            sd, sd
            % oddballs on first pitch/intensity
            sd*0.5+1/h.Settings.df*0.5, sd
            sd*1.5-1/h.Settings.df*0.5, sd
            sd*0.5+1/h.Settings.df*0.25, sd
            sd*1.5-1/h.Settings.df*0.25, sd
            % oddballs on second pitch/intensity
            sd, sd*0.5+1/h.Settings.df*0.5
            sd, sd*1.5-1/h.Settings.df*0.5
            sd, sd*0.5+1/h.Settings.df*0.25
            sd, sd*1.5-1/h.Settings.df*0.25
            ];

        %% SEQUENCE
        h.Settings.oddprob = [
            % standard
            0.8
            % oddballs on first pitch/intensity
            0.025
            0.025
            0.025
            0.025
            % oddballs on second pitch/intensity
            0.025
            0.025
            0.025
            0.025
            ]';

        h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

        % index of oddball value that are standards
        h.Settings.standardind = 1; % does not apply to 'roving oddball' design
        h.Settings.oddind = 2; % does not apply to 'roving oddball' design

        % keep oddball trials apart by at least sep_odd standards
        h.Settings.sep_odd = 2;

        % for each set, ensure a number of leading standards 
        h.Settings.std_lead = 0;

        % number of minimum sets to randomised together
        h.Settings.n_set = 1; % 1 = use min set size

        % min number of oddballs within each CP condition
        h.Settings.n_odd = [1]; % overrides h.Settings.totdur
        % min number of oddballs per randomised set, per CP
        h.Settings.n_odd_set = [1]; % overrides h.Settings.totdur

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
        %h.Settings.adaptive.startinglevel = 100; % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100)

        % number of trials of plot
        h.Settings.plottrials=0;

end

