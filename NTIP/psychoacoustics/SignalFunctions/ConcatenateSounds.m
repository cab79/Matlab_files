function soundout = ConcatenateSounds(varargin)
%
% soundout = ConcatenateSounds(s1, s2, s3, ...)
%
% This function enables to concatenate an infinite number of sounds.
%
% % EXAMPLE: create a "do, re, mi" scale
% sf = 22050;
% do = GenerateTone (sf, 261.6, 500); % this is a do4 pure tone
% re = GenerateTone (sf, 293.6, 500); % this is a re4 pure tone
% mi = GenerateTone (sf, 329.6, 500); % this is a mi4 pure tone
% doremi = ConcatenateSounds(do, re, mi);
% sound(doremi, sf)

soundout = [];
for i = 1:nargin
    soundout = [soundout; varargin{i}];
end;


