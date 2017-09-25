function [pos_ans, q] = DurationDiscriminationComplexTone_MLP(std_level, var_level, nAFC)

% Duration discrimination for a complex tone. The tone has four harmonics
% (f0=330-Hz, mi4). The subject has to tell the longest tone. The tone has
% raised cosine onset and offset gates of 10-ms.
% - PARAMETER VARIED ADAPTIVELY: the duration of the variable tone;
% - STANDARD LEVEL: here it corresponds to the duration of the standard
% tone.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
f0 = 330; % fundamental frequency of the tone in Hz
freqs = [f0*1, f0*2, f0*3, f0*4];

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateTone(sf, var_level, freqs);
variable = GenerateEnvelope(sf, variable);
standard = GenerateTone(sf, std_level, freqs);
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
q = ['Where was the longest tone (' num2str(1:nAFC) ')?: '];