function [pos_ans, q] = ProfileAnalysis_MLP(std_level, var_level, nAFC)

% Profile Analysis. In this experiment the subject listens to three complex
% tones. Two are identical (the standards). They have five harmonics all at
% the same amplitude (f0=330-Hz, mi4). The third has a similar harmonic
% structure, however, the amplitude of the third harmonic component is
% higher producing a different timbre in comparison to the standards. The
% subject has to tell the odd timbre tone. The overall level of standards
% and variable is varied randomly from trial to trial within a range of
% 5-dB. Onset and offset of tones are gated on an off with two 10-ms raised
% cosine ramps. This experiment can be run as 3AFC only. The threshold is
% given in dB.
% - PARAMETER VARIED ADAPTIVELY: the level (in dB) of the 3rd component;
% - STANDARD LEVEL: the level (in dB) of the remaining components.

if nAFC < 3
    error('This experiment must be run with three alternatives');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
f0 = 330; % fundamental frequency of the tone in Hz
dur = 250; % tones' duration in ms
freqs = [f0*1, f0*2, f0*3, f0*4, f0*5];
stand_amps = [10^(std_level/20), 10^(std_level/20), 10^(std_level/20), 10^(std_level/20), 10^(std_level/20)];
var_amps = [10^(std_level/20), 10^(std_level/20), 10^(var_level/20), 10^(std_level/20), 10^(std_level/20)];

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateTone(sf, dur, freqs, stand_amps);
standard = GenerateEnvelope(sf, standard);
standard1 = AttenuateSound(standard, -(10+rand*5));
standard2 = AttenuateSound(standard, -(10+rand*5));

variable = GenerateTone(sf, dur, freqs, var_amps);
variable = GenerateEnvelope(sf, variable);
variable = AttenuateSound(variable, -(10+rand*5));

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
isi = GenerateSilentInterval(sf, 500);
pos = randperm(3);
if pos(1) == 1
    sequence = [variable; isi; standard1; isi; standard2];
    pos_ans = 1;
elseif pos(1) == 2
    sequence = [standard1; isi; variable; isi; standard2];
    pos_ans = 2;
elseif pos(1) == 3
    sequence = [standard1; isi; standard2; isi; variable];
    pos_ans = 3;
end;

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the odd timbre tone (' num2str(1:nAFC) ')?: '];