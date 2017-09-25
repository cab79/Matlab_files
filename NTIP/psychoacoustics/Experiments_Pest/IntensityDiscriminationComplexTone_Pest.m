function [pos_ans, q] = IntensityDiscriminationComplexTone_Pest(std_level, var_level, nAFC)

% Intensity discrimination threshold for a 250-ms complex tone. The tone
% has four harmonics (f0=330-Hz, mi4). The subject has to tell the loudest
% tone. The onset and offset of the tones are gated with two 10-ms raised
% cosine ramps.
% - PARAMETER VARIED ADAPTIVELY: the level of the variable tone;
% - STANDARD LEVEL: the level of the standard tone.

if ~nAFC 
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
f0 = 330; % fundamental frequency of the tone in Hz
freqs = [f0*1, f0*2, f0*3, f0*4];
dur = 250; % tones' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateTone(sf, dur, freqs);
standard = GenerateEnvelope(sf, standard);
standard = AttenuateSound(standard, std_level);

variable = GenerateTone(sf, dur, freqs);
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
q = ['Where was the loudest tone (' num2str(1:nAFC) ')?: '];