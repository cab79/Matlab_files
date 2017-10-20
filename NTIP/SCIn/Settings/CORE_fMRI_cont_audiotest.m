function h = CORE_fMRI_cont(h,opt)

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'OddballCreateSequence';

switch opt
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'practice','fMRI_run'};
    
    case 'practice'
        
    case 'fMRI_run'
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=1;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = 1./(h.Settings.trialdur);
        % stimulus types
        % rows = conditions (stimuli randomised together), columns = outputs within each condition
        % if using D188, would be the output channels 
        % if auditory, would be the pitch, intensity or channel indices (actual channel numbers on device are
        % specified later by h.Settings.stimchan).
        h.Settings.stimtypeouts_init = [1 2; 3 4]; % [1 2; 3 4] digit 1 vs digit 5; left vs. right ear
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change condition, column 2 = change condition.
        h.Settings.cond_num_init = [1 2; 3 4]; %[1 2; 3 4]; 
        %number of output channel types, e.g. hands being stimulated, modalities (auditory plus tactile).
        h.Settings.num_chantypes = 2;
        %on each hand, number of conditions (e.g. 3 digit changes vs 1 digit changes; or two different pitch/intensity changes)
        h.Settings.num_oddballtypes = 1;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.2 0.4];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block. Sum must be integer divider of h.Settings.nchanges_per_cond
        h.Settings.cp_stims_init = [2 4 6 8; 1 2 3 4];
        %number of blocks of each change probability (CP)
        h.Settings.cond_rep_init = [2 2]; % this DIVIDES nchanges_per_cond, i.e. does not increase the number of changes
        %number of stimulus location changes within each CP condition
        h.Settings.nchanges_per_cond = [80, 80]; % multiply desired n changes PER CP CONDITION by cond_rep_init to get this number
        
        % BLOCKING/RUN OPTIONS
        % 'divide' = equally divide trials by nblocks; 
        % 'cond' = separate block for each condition
        h.Settings.blockopt = 'divide';
        % further options for 'divide':
            % number of blocks (containing multiple conditions)
            h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
            %distribute conditions equally among blocks
            h.Settings.distblocks = 1;
        % options to start sequence at beginning of every run
        % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
        % cell array
        h.Settings.blockstart = {'labjack'}; % audio,labjack,audio,
        % names of any audiofiles
        h.Settings.audiofile = {'instruct.wav','start.wav'}; % labjack
        % number of scanner triggers to wait for before starting the
        % sequence
        h.Settings.num_scanner_trig = 4;
        % wait time between scanner triggers (s)
        h.Settings.waittime_scanner_trig = 1;
        
end

function h = setgeneral(h)

%% TRIALS OR CONTINUOUS?
h.Settings.design = 'continuous';
% if continuous, how many trials ahead should be in the player schedule?
% (applies to stimulation via soundcard only)
h.Settings.ntrialsahead = 2;  %0 = all trials

%% Output options
% save sinwave from all trials as part of stim sequence file
h.Settings.savesinwave = 0;
    
%% EQUIPMENT CONTROL
%Use labjack for controlling any equipment?
h.Settings.labjack=1;
% How to control stimulator? Options: audioplayer, PsychPortAudio, labjack, spt
h.Settings.stimcontrol='PsychPortAudio';
% if using PsychPortAudio with more than 2 channels:
    %7.1 soundcard channel numbers:
    %front = 1,2 - can't use these for 7.1 without new connectors
    %C-sub = 3,4 - connect to 1/2
    %rear = 5,6 - connect to 3/4
    %side = 7,8
% total number of channels, e.g. on sound card
h.Settings.nrchannels = 2; % 8 
% Use D188
%h.Settings.D188 = 1;
% message box with OK press to start experiment?
%h.Settings.startmessage = 'Press OK to continue.';
%record EEG, NS: netstation, BV: brainvision
%h.Settings.record_EEG='NS';

%% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
% minimum total duration of stimulus sequence in seconds
%h.Settings.totdur = 60; 
% duration of trial in seconds
h.Settings.trialdur = 0;
% duration of stimulus in seconds
h.Settings.stimdur = [0.5 0.5];
% Pattern type method: intensity, pitch. Not supported: channel, duration
h.Settings.patternmethod = 'intensity';
h.Settings.patternvalue = [1 0]; % one per stimdur
% 'rand' or 'reg' spacing?
h.Settings.stimdurtype = 'reg';
% channels on stimulator to use
h.Settings.stimchan = [3 4]; % [7 8 3 4]
% Pitch/freq
h.Settings.f0 = 500; 
% Intensity
h.Settings.inten = 1; 
% Apply Tukey window?
h.Settings.Tukey = 0.25;
h.Settings.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
% sampling rate
h.Settings.fs = 96000; % don't change this, unless sure soundcard supports higher rates
h.Settings.atten = 0; % attenuation level in decibels

%% Condition-dependent stimulus parameters
% Condition method: do the stimtypes indexed by
% h.Settings.stimtypeouts_init (oddball method) differ by any other
% characteristic? E.g. different modalities may requrie different pitches
h.Settings.conditionmethod = {'pitch','intensity'};
h.Settings.conditionvalue = [200 200 500 500; 1 1 1 1];% for each number in h.Settings.stimtypeouts_init,
% Oddball method: intensity, pitch, channel
h.Settings.oddballmethod = 'channel'; % can use same type for pattern only if oddball intensity is adaptive
h.Settings.oddballvalue = [1 2, 1 2]; % [7 8 3 4]

%% SEQUENCE: For future adoption of CreateSequence
h.Settings.oddballtype = 'roving'; % options: 'roving', 'classical'
% index of oddball value that are standards
h.Settings.standardind = 1; % does not apply to 'roving oddball' design
h.Settings.oddind = 2; % does not apply to 'roving oddball' design
% keep oddball trials apart by at least sep_odd standards
h.Settings.sep_odd = 2;
% for each set, ensure a number of leading standards 
h.Settings.std_lead = 0;
% number of minimum sets to randomised together
h.Settings.n_set = 1; % 1 = use min set size

%% RESPONSE PARAMETERS
% buttonpress options: key: keyboard inputs. Blank for no button press
h.Settings.buttontype='key';
% range of keyboard presses indicating a recordable response
h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
%display RT, 1=yes, 0=no
h.Settings.displayRT=0;
% response probe
h.Settings.RPconds=[2 4]; % condition numbers to apply this
h.Settings.RPprob=0.1;
h.Settings.RPmethod='intensity';
h.Settings.RPdur = [0.15 0.2 0.15 0.5];
h.Settings.RPvalue=[1 0 1 0];

%% ADAPTIVE
% adaptive staircase: meanings of the buttonopt
%h.Settings.adaptive.buttonfun = {'low','high'}; 
% adaptive staircase: corresponding signal values that would signify a
% correct answer
%h.Settings.adaptive.signalval = [1 2];
% starting level of adaptive staircase
%h.Settings.adaptive.startinglevel = 10; % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100)

% number of trials of plot
h.Settings.plottrials=0;
