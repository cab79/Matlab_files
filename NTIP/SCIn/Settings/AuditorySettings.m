function handles = AuditorySettings(handles)

%% NOTES ON SETTINGS:
% 1/handles.df*0.25 AND (1/df+1/pitchdiff)*0.25 needs to be a multiple of 1/f0 to prevent clicking caused by phase re-setting of the f0 sound
% Settings that work:
% df=10Hz: f0=200, pitchdiff=200 
% df=25Hz: f0=200, pitchdiff=200 
% df=40Hz: f0=160, pitchdiff=160
%strategy: decide on df first, then pitchdiff (ideally a multiple of df),
%then calculate:
% 1./((1/handles.df+1/pitchdiff)*0.25 ./ [1:10])
% and choose a frequency for f0.
% seems to work if f0 and pitchdiff are both multiples of df, and f0 and
% pitchdiff are also multiples of each other?

%% ENTRAINMENT SETTINGS

% total duration of entrainment in seconds
handles.dur = 60; % set to 60 for testing, but will enentually use 600.
% sampling rate
handles.fs = 96000; % don't change this
% Left ear carrier frequency
handles.f0 = 200; 
% Binarual beats frequency: creates right ear frequency of f0+df
handles.df = 10; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.

%% ERP / FREQ TAG SETTINGS
% Note: advise only using pitch changes OR intensity changes, not both!

% Frequency of pitch changes
handles.fpitch = 2.5; % Hz - 1/fpitch must be an integer multiple of 1/handles.df. Set to 0 to turn off pitch changes.
% Pitch difference: second pitch created is f0+pitchdiff
handles.pitchdiff = 200;
% Frequency of pitch changes
handles.finten = 0; % Hz - 1/finten must be an integer multiple of 1/handles.df. Set to 0 to turn off pitch changes.
% Intensity difference: e.g 0.75 produces 75% of the intensity to
% alternate with 100% intensity.
handles.intendiff = 0.75; % multiple of normal intensity (value between 0 and 1)
% tone duration: only if distinct tones are needed with gaps of silence (not needed for current experiment)
handles.tone_dur = 0; % set to 0 for tones to fill whole trial

%% DURATION DEVIANT TRIALS SETTINGS
% i.e. all possible duration options and their probability

fc = max(handles.finten, handles.fpitch); % dont change: finds the frequency of change (of either pitch or intensity)
handles.fc=fc;
% options:
% left column = 1st inten/pitch
% middle column = 2nd inten/pitch
% right column = probability
handles.duropt = [
    % standard
    1/fc 1/fc 0.8
    % oddballs on left
    1/fc*0.5+1/handles.df*0.5, 1/fc, 0.025
    1/fc*1.5-1/handles.df*0.5, 1/fc, 0.025
    1/fc*0.5+1/handles.df*0.25, 1/fc, 0.025
    1/fc*1.5-1/handles.df*0.25, 1/fc, 0.025
    % oddballs on right
    1/fc, 1/fc*0.5+1/handles.df*0.5, 0.025
    1/fc, 1/fc*1.5-1/handles.df*0.5, 0.025
    1/fc, 1/fc*0.5+1/handles.df*0.25, 0.025
    1/fc, 1/fc*1.5-1/handles.df*0.25, 0.025
    ];

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
handles.SeqFun = 'AuditoryCreateSequence';