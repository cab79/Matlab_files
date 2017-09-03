function [SoundSequence, pos_ans] = ShuffleSounds(sf, standard, variable, nAFC, ISI, mode)
%
% [SoundSeq, pos] = ShuffleSounds(sf, standard, variable, nAFC)
% [SoundSeq, pos] = ShuffleSounds(sf, standard, variable, nAFC, ISI)
% [SoundSeq, pos] = ShuffleSounds(sf, standard, variable, nAFC, ISI, 'mode')
%
% This function enables to randomize the order of standard(s) and variable
% sounds in multiple intervals, multiple forced choice tasks. Sounds are
% concatenated into a single sound vector (SoundSeq) and separated by
% silent intervals of ISI duration. The function returns the position (pos)
% of the variable sound. ShuffleSounds works with either monophonic and
% stereophonic sounds.
% 
% SF: sample frequency of the sounds (in Hz)
% STANDARD: sound vector containing the "standard" sound
% VARIABLE: sound vector containing the "variable" sound
% nAFC: number of altervative foced choice
% ISI: duration (in ms) of the silent interval between the sounds (default
% is 500-ms)
% MODE: a string. Permitted values are 'AXB', 'XAB' or 'random' (default).
% It can be used in 2AFC tasks only. In these tasks the user can decide to
% present a sequence of three sounds (two standards and one variable) with
% one of the two standard sounds (i.e., X) being always the middle sound of
% the sound triplet (AXB) or the first sound of the sound triplet (XAB).

if nargin == 5, mode='random'; end;
if nargin == 4, mode='random'; ISI = 500; end;
if nargin == 3, mode='random'; ISI = 500; nAFC = 2; end;

SoundSequence = [];
% generate the ISI silent interval
if size(standard, 2)==2
    isi_silence = zeros(round(sf*(ISI/1000)), 2);
    shortsilence = zeros(100, 2);
else
    isi_silence = zeros(round(sf*(ISI/1000)), 1);
    shortsilence = zeros(100, 1);
end
    
switch mode
    case 'random'
        soundsposition = randperm(nAFC);
        for i=1:nAFC
            if soundsposition(i)==1
                SoundSequence = [SoundSequence; variable];
                pos_ans = i;
            else
                SoundSequence = [SoundSequence; standard];
            end
            if i<nAFC
                SoundSequence = [SoundSequence; isi_silence];
            end
        end
    case 'AXB'
        if nAFC ~=2, error('nAFC must be equal to "2"'); end;
        soundsposition = randperm(nAFC);
        for i=1:nAFC
            if soundsposition(i)==1
                SoundSequence = [variable; isi_silence; standard; isi_silence; standard];
                pos_ans = 1;
            else
                SoundSequence = [standard; isi_silence; standard; isi_silence; variable];
                pos_ans = 2;
            end
        end
    case 'XAB'
        if nAFC ~=2, error('nAFC must be equal to "2"'); end;
        soundsposition = randperm(nAFC);
        for i=1:nAFC
            if soundsposition(i)==1
                SoundSequence = [standard; isi_silence; variable; isi_silence; standard];
                pos_ans = 1;
            else
                SoundSequence = [standard; isi_silence; standard; isi_silence; variable];
                pos_ans = 2;
            end
        end
end

% add a 100 samples silence before and after the SoundSequence to prevent
% clicks in the  on/off switching of the sound card.
SoundSequence =[shortsilence; SoundSequence; shortsilence];