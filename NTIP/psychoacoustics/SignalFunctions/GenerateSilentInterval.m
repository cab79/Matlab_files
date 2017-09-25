function silence = GenerateSilentInterval(sf, dur, stereo);
%
% silence = GenerateSilentInterval(sf, dur)
% silence = GenerateSilentInterval(sf, dur, 'stereo')
%
% This function generates a silent interval or the desired duration.
% Default output is a monophonic silent interval.
%
% SF: sample frequency of the silence in Hz
% DUR: duration of the silent interval in ms
% STEREO: type stereo as a string (i.e., 'stereo') if you need a
% stereophonic silence

if nargin == 3
    silence = zeros(round(sf*(dur/1000)), 2);
elseif nargin <3
    silence = zeros(round(sf*(dur/1000)), 1);
end