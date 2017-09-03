function [pos_ans, q] = MelodyMistuningDetection_MLP(std_level, var_level, nAFC)

% Melody mistuning detection. The major diatonic equitempered scale (i.e.,
% do, re mi, fa, sol, la, ti, do) is played (starting do, do4=261.6-Hz).
% The sol note has a variable pitch. The subject has to tell whether the
% scale is in tune or out of tune (in yes/no task) or to tell the out of
% tune scale (in nAFC task). Notes are 500-ms complex tones of five
% harmonics. All tones are gated on and off with two raised cosine ramps of
% 10-ms. The threshold is estimated in cents. To convert the threshold in
% hertz: threshold=261.6*2^((700+t)/1200). Where t is the estimated
% threshold in cents.
% - PARAMETER VARIED ADAPTIVELY: the frequency in cents of the 5th note of
% the scale;
% - STANDARD LEVEL: here it corresponds to the catch trial level.

%%% BEGINNING EXPERIMENT'S PARAMETER %%%
sf = 44100; % sample frequency in Hz
dur = 400; % tones' duration in ms
f_do_start = 261.6; % fundamental frequency of the tone in Hz and first note of the scale
f_re = f_do_start * 2^(200/1200);
f_mi = f_do_start * 2^(400/1200);
f_fa = f_do_start * 2^(500/1200);
f_sol = f_do_start * 2^((700+std_level)/1200);
f_la = f_do_start * 2^(900/1200);
f_ti = f_do_start * 2^(1100/1200);
f_do_end = f_do_start * 2^(1200/1200);


% [1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SOUNDS
do1 = GenerateTone(sf, dur, [f_do_start*1, f_do_start*2, f_do_start*3, f_do_start*4]);
re = GenerateTone(sf, dur, [f_re*1, f_re*2, f_re*3, f_re*4]);
mi = GenerateTone(sf, dur, [f_mi*1, f_mi*2, f_mi*3, f_mi*4]);
fa = GenerateTone(sf, dur, [f_fa*1, f_fa*2, f_fa*3, f_fa*4]);
sol = GenerateTone(sf, dur, [f_sol*1, f_sol*2, f_sol*3, f_sol*4]);
la = GenerateTone(sf, dur, [f_la*1, f_la*2, f_la*3, f_la*4]);
ti = GenerateTone(sf, dur, [f_ti*1, f_ti*2, f_ti*3, f_ti*4]);
do2 = GenerateTone(sf, dur, [f_do_end*1, f_do_end*2, f_do_end*3, f_do_end*4]);

% the variable parameter
f_sol_var = f_do_start * 2^((700+var_level)/1200);
sol_var = GenerateTone(sf, dur, [f_sol_var*1, f_sol_var*2, f_sol_var*3, f_sol_var*4]);

do1 = GenerateEnvelope(sf, do1);
re = GenerateEnvelope(sf, re);
mi = GenerateEnvelope(sf, mi);
fa = GenerateEnvelope(sf, fa);
sol = GenerateEnvelope(sf, sol);
la = GenerateEnvelope(sf, la);
ti = GenerateEnvelope(sf, ti);
do2 = GenerateEnvelope(sf, do2);

sol_var = GenerateEnvelope(sf, sol_var);

variable = ConcatenateSounds(do1, re, mi, fa, sol_var, la, ti, do2);
standard = ConcatenateSounds(do1, re, mi, fa, sol, la, ti, do2);

% [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOMIZE POSITION OF STANDARD AND VARIABLE AND SET, ACCORDINGLY, THE KEY
% THE SUBJECT HAS TO PRESS TO GIVE A CORRECT RESPONSE
if nAFC == 0
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
    q = ['Is the scale out of tune ("1") or in tune ("0")?: '];
else
    q = ['In which position is the mistuned scale? (' num2str(1:nAFC) ')?: '];
end;