function tone = GenerateSweep(sf, dur, f_start, f_end, mode, amps, phases)
%
% tone = GenerateSweep(sf, dur, f_start, f_end)
% tone = GenerateSweep(sf, dur, f_start, f_end, mode)
% tone = GenerateSweep(sf, dur, f_start, f_end, mode, amps)
% tone = GenerateSweep(sf, dur, f_start, f_end, mode, amps, phases)
%
% This function generates a frequency-sweep tone starting from a frequency
% f_start and ending at a frequency f_end. The tone can be either pure or
% complex, the sweep can be either exponential (default) or linear. In the
% case of complex frequency-sweeps the user can also specify the amplitude
% of each harmonic of the complex tone.
%
% SF: sample frequency of the tone in Hz
% DUR: duration of the tone in ms
% F_START: the starting frequency (or the frequencies) of the sweep.
% F_END: the ending frequency (or the frequencies) of the sweep.
% MODE: the sweep type, either 'exponential' or 'linear'
% AMPS: a variable containing the linear amplitude of each frequency
% component (default is 1)
% PHASES: a variable containing the phases of the frequency components in
% radians (default is 0)
%
% % EXAMPLE: generate a pure tone frequency sweep starting from 5000-Hz and
% descending exponentially in frequency to 100-Hz in 3000-ms
% tone = GenerateSweep(44100, 3000, 5000, 100);

if nargin == 6, phases=zeros(size(f_start)); end
if nargin == 5, phases=zeros(size(f_start)); amps=ones(size(f_start)); end;
if nargin == 4, phases=zeros(size(f_start)); amps=ones(size(f_start)); mode='exponential'; end;

numberofsamples = sf * dur/1000;
t = (0:numberofsamples-1) / sf;

switch mode
    case 'exponential'
        tone = zeros(size(t));
        % calculate che exponent of the sweep
        k = zeros(size(f_start));
        k = (f_end./f_start).^(1/(dur/1000));

        % generate and add each frequency component
        for i=1:length(f_start)
            f_component = amps(i) * sin(phases(i) + 2*pi*f_start(i).*(k(i).^t - 1)/log(k(i)));
            tone = tone + f_component;
        end
        % amplitude normalisation
        tone = tone' / max(abs(tone));
        % little correction
        tone = tone * .999;
    case 'linear'
        tone = zeros(size(t));
        % calculate che exponent of the sweep
        k = zeros(size(f_start));
        k = (f_end-f_start)/(dur/1000);

        % generate and add each frequency component
        for i=1:length(f_start)
            f_component = amps(i) * sin(phases(i) + 2*pi*(f_start.*t + (k/2).*t.^2));
            tone = tone + f_component;
        end
        % amplitude normalisation
        tone = tone' / max(abs(tone));
        % little correction
        tone = tone * .999;
end