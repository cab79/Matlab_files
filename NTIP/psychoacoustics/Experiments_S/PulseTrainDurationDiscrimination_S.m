function [pos_ans, q] = PulseTrainDurationDiscrimination_S(std_level, var_level, nAFC)

% Pulse-train discrimination. The standard stimulus consists of six 20-ms
% pulses of a 1-kHz tone. These pulses are arranged in three pairs, with 40
% ms of silence between members of a pair and 120 ms between pairs. The
% temporal structure of the variable sequence is varied by increasing the
% separation between members of each pair, with a corresponding decrease in
% the between-pair time and, thus, a constant interval between the first
% tones in each of the successive pairs. Thus, the first, third, and fifth
% tones are fixed in time, while the onsets of the second, fourth, and
% sixth tones are delayed by varying amounts.
% - PARAMETER VARIED ADAPTIVELY: the duration of the gap within/between
% pairs of the variable stimulus;
% - STANDARD LEVEL: in the current experiment it is not in use

if nAFC < 3
    error('This experiment requires a minimum of three alternatives');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
pulse_duration = 20;
interval_between_pulses_s = 120;
interval_within_pulses_s = 40;
interval_between_pulses_v = 120-var_level;
interval_within_pulses_v = 40+var_level;
pulse_freq = 1000;

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
pulse = GenerateTone(sf, pulse_duration, pulse_freq);
pulse = GenerateEnvelope(sf, pulse);
ss1 = GenerateSilentInterval(sf, interval_between_pulses_s);
ss2 = GenerateSilentInterval(sf, interval_within_pulses_s);
sv1 = GenerateSilentInterval(sf, interval_between_pulses_v);
sv2 = GenerateSilentInterval(sf, interval_within_pulses_v);
standard = ConcatenateSounds(pulse, ss2, pulse, ss1, pulse, ss2, pulse, ss1, pulse, ss2, pulse);
variable = ConcatenateSounds(pulse, sv2, pulse, sv1, pulse, sv2, pulse, sv1, pulse, sv2, pulse);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the odd rhythm (' num2str(1:nAFC) ')?: '];