function [pos_ans, q] = PitchDiscriminationComplexTone_S(std_level, var_level, nAFC)

% Pitch discrimination threshold for a 250-ms complex tone. The tone has
% four harmonics. The subject has to tell the highest pitch tone. Onset and
% offset of tones are gated on and off with two 10-ms raised cosine ramps.
% - PARAMETER VARIED ADAPTIVELY: the delta: the frequency difference
% between the standard and the variable tone;
% - STANDARD LEVEL: here it corresponds to the frequency of standard tone.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
standard_freqs = [std_level*1, std_level*2, std_level*3, std_level*4];
dur = 250; % tones' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateTone(sf, dur, [(std_level+var_level)*1, (std_level+var_level)*2, (std_level+var_level)*3, (std_level+var_level)*4]);
variable = GenerateEnvelope(sf, variable);
standard = GenerateTone(sf, dur, standard_freqs);
standard = GenerateEnvelope(sf, standard);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the highest pitch tone (' num2str(1:nAFC) ')?: '];