function h = ASOT(h,opt)

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'OddballCreateSequence';

switch opt
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'1'};
    
    case '1'
        % set general options
        h = setgeneral(h);
        %record response, 1=yes, 0=no
        h.Settings.record_response=1;
        %between-train frequency (1000 divided by ISI in ms)
        h.Settings.freq = 1./(h.Settings.trialdur);
        % stimulus types
        % rows = conditions (stimuli randomised together), columns = output points within each condition
        % if using D188, would be the output channels 
        % if auditory, would be the pitch or intensity options
        h.Settings.stimtypeouts_init = [1 2];
        %number of hands being stimulated
        h.Settings.num_hands = 1;
        %condition numbers to associate with each condition; rows = conditions,
        %column1 = no change conditions, column 2 = change condition.
        h.Settings.cond_num_init = [1 2];
        %on each hand, number of conditions (e.g. 3 digit changes vs 1 digit changes; or two different pitch/intensity changes)
        h.Settings.num_changes = 1;
        %list of probabilities of a stimulus changing location on each trial
        h.Settings.change_prob_init = [0.2 0.4];
        %for each change probability (rows), list of number of trials in which the stimuli remain in same
        %location - these will be randomised within a block. Sum must be integer divider of h.Settings.nchanges_per_cond
        h.Settings.cp_stims_init = [2 4 6 8; 1 2 3 4];
        %number of blocks of each change probability
        h.Settings.cond_rep_init = [1 1];
        %number of stimulus location changes within each condition
        h.Settings.nchanges_per_cond = [20, 20]; 
        
        % BLOCKING OPTIONS
        % 'divide' = equally divide trial by nblocks; 
        % 'cond' = separate block for each condition
        h.Settings.blockopt = 'cond';
        % further options for 'divide':
            % number of blocks (containing multiple conditions) in which experiment is run until a pause
            h.Settings.nblocks = 2; 
            %distribute conditions equally among blocks
            h.Settings.distblocks = 0;
        % options to start sequence at beginning of every block
        % 'msgbox', 'labjack', 'buttonpress'
        h.Settings.blockstart = 'buttonpress';
        % require additional button press to start and re-start, including after a pause? 0 or 1
        h.Settings.buttonstart = 0;
        
end

function h = setgeneral(h)

%% TRIALS OR CONTINUOUS?
h.Settings.design = 'trials';

%% EQUIPMENT CONTROL
%Use labjack for controlling any equipment?
%h.Settings.labjack=1;
% How to control stimulator? Options: audioplayer, PsychPortAudio, labjack, spt
h.Settings.stimcontrol='PsychPortAudio';
% if using PsychPortAudio with more than 2 channels:
    %7.1 soundcard channel numbers:
    %front = 1,2
    %C-sub = 3,4
    %rear = 5,6
    %side = 7,8
h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
h.Settings.stimchan = [1 2]; % channels on stimulator to use
% Use D188
%h.Settings.D188 = 1;
% message box with OK press to start experiment?
%h.Settings.startmessage = 'Press OK to continue.';
%record EEG, NS: netstation, BV: brainvision
%h.Settings.record_EEG='NS';

%% STIMULUS PARAMETERS
% Oddball method: intensity, pitch, channel
h.Settings.oddball = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
% Pattern type method: intensity, pitch, channel, duration
h.Settings.pattern = 'pitch';
% duration of stimulus sequence in seconds
%h.Settings.totdur = 60; 
% duration of trial in seconds
h.Settings.trialdur = 0.5;
% duration of stimulus in seconds
%h.Settings.stimdur = 0.25; 
h.Settings.stimdur = [0.02:0.02:0.10]; % could also randomise these
% 'rand' or 'reg' spacing?
h.Settings.stimdurtype = 'rand';
%REMOVE: h.Settings.npulses = 1;
% Apply Tukey window?
h.Settings.Tukey = 0.25;
% ratio of on:off periods during stimulu
h.Settings.ratio_on = 1; % not yet implemented
% sampling rate
h.Settings.fs = 96000; % don't change this, unless sure soundcard supports higher rates
% Pitch
h.Settings.f0 = [500:100:900]; % uses only first, or both if oddball selected
% Intensities
h.Settings.inten = [1]; % one per column of h.Settings.stimdur OR oddball condition
h.Settings.atten = -30; % attenuation level in decibels
%% RESPONSE PARAMETERS
% buttonpress options: key: keyboard inputs. Blank for no button press
h.Settings.buttontype='key';
% range of keyboard presses indicating a recordable response
h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
%display RT, 1=yes, 0=no
h.Settings.displayRT=0;

%% ADAPTIVE
% adaptive staircase: meanings of the buttonopt
h.Settings.adaptive.buttonfun = {'low','high'}; 
% adaptive staircase: corresponding signal values that would signify a
% correct answer
h.Settings.adaptive.signalval = [1 2];
% starting level of adaptive staircase
h.Settings.adaptive.startinglevel = 10; % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100)

