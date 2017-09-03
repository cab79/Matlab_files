function [pos_ans, q] = GapDetectionWhiteNoise_S(std_level, var_level, nAFC)

% Gap detection. A band of 750-ms gaussian noise has a gap in its temporal
% center. Gap duration is varied according to the listener performance. The
% noise has 0.5-ms cosine ramps at the beginning and end of the gap. In
% nI-nAFC tasks, the standard is always a 750-ms broadband noise with no
% gap whereas the variable contains the gap.
% - PARAMETER VARIED ADAPTIVELY: the duration of the gap;
% - STANDARD LEVEL: here it corresponds to the overall noise duration.

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateNoise(sf, std_level);
% split the standard in two portions
noise1 = standard(1:round((end/2)), 1);
noise2 = standard(round((end/2)+1):end, 1);
gapdurationinsamples = floor((var_level/1000) * sf);
noise1 = noise1(1:(end-round(gapdurationinsamples/2)), 1);
noise2 = noise2(round(gapdurationinsamples/2):end, 1);
gap = GenerateSilentInterval(sf, var_level);
standard = GenerateEnvelope(sf, standard);
noise1 = GenerateEnvelope(sf, noise1, 10, 0.5);
noise2 = GenerateEnvelope(sf, noise2, 0.5, 10);
variable = ConcatenateSounds(noise1, gap, noise2);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
if ~nAFC 
    pos_ans = 1;
    sequence = variable;
else
    [sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);
end;

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
if nAFC == 0
    q = 'Can you hear ("1") the gap or not ("0")?: ';
else
    q = ['Where was the noise with the gap (' num2str(1:nAFC) ')?: '];
end;