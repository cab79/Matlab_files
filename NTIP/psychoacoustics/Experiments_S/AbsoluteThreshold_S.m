function [pos_ans, q] = AbsoluteThreshold_S(std_level, var_level, nAFC)

% Absolute threshold for a 500-ms pure tone of 1-kHz. The tone is gated on
% and off with two raised cosine ramps of 10-ms.
% - PARAMETER VARIED ADAPTIVELY: the level of the tone;
% - STANDARD LEVEL: in the current experiment this parameter has no use;
% !!!PLEASE NOTE THAT THE EXACT SOUND PRESSURE LEVEL AT LISTENER'S EAR MUST
% BE DETERMINED APART BY MEANS OF AN ADEQUATE HARDWARE!!!

if nAFC
    error('This experiment is for "yes/no" task only!');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 500; % overall duration of noise
freq = 1000; % tone frequency

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
variable = GenerateTone(sf, dur, freq);
variable = GenerateEnvelope(sf, variable);
variable = AttenuateSound(variable, var_level);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
pos_ans = 1;

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(variable, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = 'Can you hear the tone ("1") or not ("0")?: ';