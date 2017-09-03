function [pos_ans, q] = IntensityDiscriminationWhiteNoise_MLP(std_level, var_level, nAFC)

% Intensity discrimination threshold for a 250-ms white noise. The subject
% has to tell the loudest noise. The onset and offset of the noises are
% gated with two 10-ms raised cosine ramps.
% - PARAMETER VARIED ADAPTIVELY: the level of the variable noise;
% - STANDARD LEVEL: here it corresponds to the level of the standard noise.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

% EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 250; % noise' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateNoise(sf, dur);
standard = GenerateEnvelope(sf, standard);
standard = AttenuateSound(standard, std_level);

variable = GenerateNoise(sf, dur);
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
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the loudest noise (' num2str(1:nAFC) ')?: '];