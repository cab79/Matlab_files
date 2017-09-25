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

%if ~exist('opt','var')
%    opt = {'10Hz'};
%end

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'CreateSequence';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'10Hz'};
    
    case '10Hz'
        
    %% TRIALS OR CONTINUOUS?
    h.Settings.design = 'continuous';
    
    %% ENTRAINMENT SETTINGS

    % duration of trial in seconds
    h.Settings.trialdur = 0; % 0 = all stimuli are continuous
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 60; 
    % sampling rate
    h.Settings.fs = 96000; % don't change this
    % Left ear carrier frequency (pitch)
    h.Settings.f0 = [200 400]; % one per column of h.Settings.stimdur OR oddball condition
    % Binarual beats frequency: creates right ear frequency of f0+df
    h.Settings.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.

    %% ERP / FREQ TAG SETTINGS
    % Note: advise only using pitch changes OR intensity changes, not both!

    % Oddball type: intensity, pitch, channel, duration
    h.Settings.oddball = 'duration';
    % Pattern type method: intensity, pitch, channel, duration
    h.Settings.pattern = 'pitch';
    % Frequency of pitch changes
    h.Settings.fpattern = 2.5; % Hz - 1/fpattern must be an integer multiple of 1/h.Settings.df. Set to 0 to turn off pitch changes.
    % Intensity difference: e.g 0.75 produces 75% of the intensity to
    % alternate with 100% intensity.
    h.Settings.inten = [1 0.1]; % multiple of normal intensity (value between 0 and 1)
    h.Settings.atten = -30; % attenuation level in decibels
    % tone duration: only if distinct tones are needed with gaps of silence (not needed for current experiment)
    %h.Settings.tone_space = 0; % set to 0 for tones to fill whole trial
    

    %% DURATION DEVIANT TRIALS SETTINGS
    % i.e. all possible duration options and their probability

    fc = h.Settings.fpattern; % dont change: finds the frequency of change (of either pitch or intensity)
    %h.Settings.fc=fc;
    
    %STIMDUR
    % In general, each row is a different stim type, with columns providing
    % within-stimulus temporal patterns of pitch, intensity or channel.
    % Here the options are chosen as:
        % left column = 1st inten/pitch
        % right column = 2nd inten/pitch
    % and the temporal pattern is defined by fc (from either fpitch or finten)
    h.Settings.stimdur = [
        % standard
        1/fc, 1/fc
        % oddballs on first pitch/intensity
        1/fc*0.5+1/h.Settings.df*0.5, 1/fc
        1/fc*1.5-1/h.Settings.df*0.5, 1/fc
        1/fc*0.5+1/h.Settings.df*0.25, 1/fc
        1/fc*1.5-1/h.Settings.df*0.25, 1/fc
        % oddballs on second pitch/intensity
        1/fc, 1/fc*0.5+1/h.Settings.df*0.5
        1/fc, 1/fc*1.5-1/h.Settings.df*0.5
        1/fc, 1/fc*0.5+1/h.Settings.df*0.25
        1/fc, 1/fc*1.5-1/h.Settings.df*0.25
        ];
    
    % 'rand' or 'reg' spacing?
    h.Settings.stimdurtype = 'reg'; % not needed unless 'rand'
    
    h.Settings.durprob = [
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
        ];
    
    %% Equipment options
    % record EEG, NS: netstation, BV: brainvision
    h.Settings.record_EEG='';
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'left','right'}; 
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % require button press to start and re-start after a pause? 0 or 1
    h.Settings.buttonstart = 1;
    % Use labjack for controlling any equipment?
    h.Settings.labjack=0;
    % How to control stimulator? Options: audioplayer, labjack, spt
    h.Settings.stimcontrol='audioplayer';
    h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stimchan = [1 2]; % channels on stimulator to use
    
    % save sinwave from all trials as part of stim sequence file
    h.Settings.savesinwave = 1;
    

end

