function [pos_ans, q] = ITD_DiscriminationPureTone_Pest(std_level, var_level, nAFC)

% ITD discrimination for a 330-Hz, 250-ms pure tone. The subject has to
% tell whether the leading tone of the tone pair comes from the left or
% from the right. Both tones have a certain ITD. This task is identical to
% Saberi and Green (P&P, 1997).
% - PARAMETER VARIED ADAPTIVELY: the ITD of the variable tone (here it is
% the left tone). The same ITD value is also used for the standard tone
% however, with opposite sign;
% - STANDARD LEVEL: in the current experiment it has no use.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
freq = 330; % frequency of the standard tone in Hz
amp = 1; % amplitude of the sine wave
phase = 0; % phase of the sine wave
gatedur = 2.5; % onset & offset gates duration  in ms
dur = 250; % tones' duration in ms

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateTone(sf, dur, freq, amp, phase, -var_level);
standard = GenerateEnvelope(sf, standard, gatedur);
standard = AttenuateSound(standard, -10);

variable = GenerateTone(sf, dur, freq, amp, phase, var_level);
variable = GenerateEnvelope(sf, variable, gatedur);
variable = AttenuateSound(variable, -10);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECT
q = ['Which was the tone coming from the left (' num2str(1:nAFC) ')?: '];