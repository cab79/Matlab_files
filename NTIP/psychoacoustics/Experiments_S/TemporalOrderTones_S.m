function [pos_ans, q] = TemporalOrderTones_S(std_level, var_level, nAFC)

% Temporal order for tones. The task is to discriminate the order in which
% two equal-duration pure tones of 550 and 710-Hz are presented. The
% duration of the two tones is varied according to listener performance.
% Tones are presented without a gap between them and are preceded and
% followed, without gaps, by 100-ms ?leader? and ?trailer? tones at 625-Hz.
% The onset and offset of the tones are gated with two 10-ms raised cosine
% ramps.
% - PARAMETER VARIED ADAPTIVELY: the duration of the middle tones;
% - STANDARD LEVEL: in the current experiment it is not in use.

if nAFC < 3
    error('We recomment to estimate this threshold with a 3AFC tasks');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
freq = 625; % frequency of leader and trailer
f1 = 550; % frequency of one tone
f2 = 710; % frequency of the other
dur = 100; % duration of leader and trailer

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
leadtrail = GenerateTone(sf, dur, freq);
leadtrail = GenerateEnvelope(sf, leadtrail);

tone1 = GenerateTone(sf, var_level, f1);
tone1 = GenerateEnvelope(sf, tone1);
tone2 = GenerateTone(sf, var_level, f2);
tone2 = GenerateEnvelope(sf, tone2);

pos = randperm(2);
if pos(1)==1
    standard = ConcatenateSounds(leadtrail, tone1, tone2, leadtrail);
    variable = ConcatenateSounds(leadtrail, tone2, tone1, leadtrail);
else
    standard = ConcatenateSounds(leadtrail, tone2, tone1, leadtrail);
    variable = ConcatenateSounds(leadtrail, tone1, tone2, leadtrail);
end;

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Which tone sequence was the odd one (' num2str(1:nAFC) ')?: '];