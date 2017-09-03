function [pos_ans, q] = DurationDiscriminationPureTone_Pest(std_level, var_level, nAFC)

% Duration discrimination for a 1-kHz pure tone. The subject has to
% tell the longest tone. The tone has raised cosine onset and offset gates
% of 10-ms.
% - PARAMETER VARIED ADAPTIVELY: the duration of the variable tone;
% - STANDARD LEVEL: the duration of the standard tone.

if ~nAFC 
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
freq = 1000; % frequency of the standard tone in Hz

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateTone(sf, var_level, freq);
variable = GenerateEnvelope(sf, variable);
standard = GenerateTone(sf, std_level, freq);
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
q = ['Where was the longest tone? (' num2str(1:nAFC) ')?: '];