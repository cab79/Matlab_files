function [pos_ans, q] = IntensityDiscriminationPureTone_Pest(std_level, var_level, nAFC)

% Intensity discrimination threshold for a 1-kHz, 250-ms pure tone. The
% subject has to tell the loudest tone. The onset and offset of the tones
% are gated with two 10-ms raised cosine ramps.
% - PARAMETER VARIED ADAPTIVELY: the level of the variable tone;
% - STANDARD LEVEL: the level of the standard tone.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
freq = 1000; % frequency of the standard tone in Hz
dur = 250; % tones' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateTone(sf, dur, freq);
standard = GenerateEnvelope(sf, standard);
standard = AttenuateSound(standard, std_level);

variable = GenerateTone(sf, dur, freq);
variable = GenerateEnvelope(sf, variable);
variable = AttenuateSound(variable, var_level);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECT
q = ['Where was the loudest tone (' num2str(1:nAFC) ')?: '];
    