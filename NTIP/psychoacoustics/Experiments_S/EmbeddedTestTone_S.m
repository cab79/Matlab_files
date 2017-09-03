function [pos_ans, q] = EmbeddedTestTone_S(std_level, var_level, nAFC)

% Subjects listen for one member of a sequence of nine tones with
% frequencies ranging from 300 to 3000-Hz. A different, randomly selected
% series of nine tones is presented on each trial. The task is to detect
% the presence of the fifth tone in the sequence. The tone is absent in the
% standard. The duration of all tones except the fifth, or target tone, is
% 40-ms. All tones have 2.5-ms raised cosine onset and offset gates. The
% test is made more difficult by reducing the duration of the target tone.
% - PARAMETER VARIED ADAPTIVELY: the duration of the middle tone;
% - STANDARD LEVEL: the duration of the remaining tones.

if nAFC<3
    error('We suggest to run this experiment with at least three alternatives');
end;

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
gatedur = 2.5; % duration of the onset and offset gates of all sounds

% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
for i=1:9
    f(i)=(rand*2700)+300; %generate the random frequencies for the tones
end;

t1 = GenerateTone(sf, std_level, f(1));
t2 = GenerateTone(sf, std_level, f(2));
t3 = GenerateTone(sf, std_level, f(3));
t4 = GenerateTone(sf, std_level, f(4));
t5 = GenerateTone(sf, var_level, f(5));
t6 = GenerateTone(sf, std_level, f(6));
t7 = GenerateTone(sf, std_level, f(7));
t8 = GenerateTone(sf, std_level, f(8));
t9 = GenerateTone(sf, std_level, f(9));

t1 = GenerateEnvelope(sf, t1, gatedur);
t2 = GenerateEnvelope(sf, t2, gatedur);
t3 = GenerateEnvelope(sf, t3, gatedur);
t4 = GenerateEnvelope(sf, t4, gatedur);
t5 = GenerateEnvelope(sf, t5, gatedur);
t6 = GenerateEnvelope(sf, t6, gatedur);
t7 = GenerateEnvelope(sf, t7, gatedur);
t8 = GenerateEnvelope(sf, t8, gatedur);
t9 = GenerateEnvelope(sf, t9, gatedur);

standard = ConcatenateSounds(t1, t2, t3, t4, t6, t7, t8, t9);
variable = ConcatenateSounds(t1, t2, t3, t4, t5, t6, t7, t8, t9);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
[sequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC);

% [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAY THE SOUND
sound(sequence, sf, 16);

% [4] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASK THE QUESTION TO THE SUBJECTs
q = ['Where was the odd tone sequence?  (' num2str(1:nAFC) ')?: '];