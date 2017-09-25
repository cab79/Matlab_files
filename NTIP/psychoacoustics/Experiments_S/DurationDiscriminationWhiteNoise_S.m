function [pos_ans, q] = DurationDiscriminationWhiteNoise_S(std_level, var_level, nAFC)

% Duration discrimination for a white noise. The subjects has to tell the
% longest noise. The noise has raised cosine onset and offset gates of
% 10-ms.
% - PARAMETER VARIED ADAPTIVELY: the delta: the difference in duration
% between variable and standard noise;
% - STANDARD LEVEL: here it corresponds to the duration of the standard
% noise.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

% EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateNoise(sf, std_level+var_level);
variable = GenerateEnvelope(sf, variable);
standard = GenerateNoise(sf, std_level);
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
q = ['Where was the longest noise?: (' num2str(1:nAFC) ')?: '];