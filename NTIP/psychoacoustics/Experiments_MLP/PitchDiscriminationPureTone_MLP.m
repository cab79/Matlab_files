function [pos_ans, q] = PitchDiscriminationPureTone_MLP(std_level, var_level, nAFC)

% Pitch discrimination threshold for a 250-ms long pure tone. The highest
% pitch tone has to be discriminated among a number of choices. Onset and
% offset of tones are gated on and off with two 10-ms raised cosine ramps.
% - PARAMETER VARIED ADAPTIVELY: the frequency of the variable tone;
% - STANDARD LEVEL: here it corresponds to the frequency of standard tones.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 250; % tones' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateTone(sf, dur, var_level);
variable = GenerateEnvelope(sf, variable);
standard = GenerateTone(sf, dur, std_level);
standard = GenerateEnvelope(sf, standard);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECT
q = ['Discriminate the highest pitch tone''s temporal position (' num2str(1:nAFC) ') '];