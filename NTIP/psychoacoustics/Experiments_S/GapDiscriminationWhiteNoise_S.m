function [pos_ans, q] = GapDiscriminationWhiteNoise_S(std_level, var_level, nAFC)

% Gap-duration discrimination. The standard is a 750-ms Gaussian noise with
% a silent of standard duration placed at its temporal center. The variable
% has a variable gap duration and the length of its gap is changed as
% function of the subject performance. All noises have a 0.5-ms cosine ramp
% at onset and offset.
% - PARAMETER VARIED ADAPTIVELY: the delta: the difference between the
% duration of gap in the standard tone and the duration of the gap in the
% variable tone;
% - STANDARD LEVEL: here it corresponds to the duration of the
% standard-noise gap.

if ~nAFC
    error('Discrimination thresholds must be estimated with nAFC tasks!');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 750/2; % overall duration of noise
isi = 800; % duration of the silent interval between the stimuli

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
noise1 = GenerateNoise(sf, dur);
noise2 = GenerateNoise(sf, dur);
gap = GenerateSilentInterval(sf, std_level);
noise1 = GenerateEnvelope(sf, noise1, 0.5, 0.5);
noise2 = GenerateEnvelope(sf, noise2, 0.5, 0.5);
standard = ConcatenateSounds(noise1, gap, noise2);

vargapdur = std_level + var_level;
vargap = GenerateSilentInterval(sf, vargapdur);
variable = ConcatenateSounds(noise1, vargap, noise2);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC, isi);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the noise with the longest gap (' num2str(1:nAFC) ')?: '];