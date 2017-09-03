function [pos_ans, q] = SAM_Detection_8Hz_S(std_level, var_level, nAFC)

% Sinusoidal Amplitude Modulation (SAM) noise discrimination. A 500-ms
% Gaussian noise is sinusoidally amplitude modulated at 8-Hz. The depth
% of the modulation is expressed as 20log(m), where m is a modulation index
% that ranges from 0.0 (no modulation) to 1.0 (full modulation). The
% subject has to detect the modulation (in yes/no task) or to tell which
% interval has the modulated noise. Modulated and unmodulated stimuli are
% equated for total RMS power. Noises have two 10-ms raised cosine ramps at
% onset and offset. The threshold is the modulation depth (in dB).
% - PARAMETER VARIED ADAPTIVELY: the modulation depth of the variable noise
% expressed as 20*log(m) where m is the modulation index that ranges from 0
% (no mudulation) to 1 (full modulation);
% - STANDARD LEVEL: in the current experiment it has no use.

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 500; % duration of the noise
mrate = 8; % modulation rate

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
standard = GenerateNoise(sf, dur);
standard = ApplyAmplitudeModulation(sf, standard, mrate, exp(std_level/20));
variable = ApplyAmplitudeModulation(sf, standard, mrate, exp(var_level/20));
standard = AttenuateSound(standard, -20);
variable = AttenuateSound(variable, -20);
standard = GenerateEnvelope(sf, standard);
variable = GenerateEnvelope(sf, variable);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
if ~nAFC
    pos_var = 1;
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
    q = 'Is the noise sinusoidally modulated in amplitude ("1") or ("0") not?: ';
else
    q = ['Where is the modulated noise (' num2str(1:nAFC) ')?: '];
end;