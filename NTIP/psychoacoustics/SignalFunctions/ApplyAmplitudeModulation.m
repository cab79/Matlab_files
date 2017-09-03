function soundout = ApplyAmplitudeModulation(sf, soundin, mrate, mdepth, mphase, ITD);
%
% soundout = ApplyAmplitudeModulation(sf, soundin, mrate, mdepth);
% soundout = ApplyAmplitudeModulation(sf, soundin, mrate, mdepth, mphase);
% soundout = ApplyAmplitudeModulation(sf, soundin, mrate, mdepth, mphase, ITD);
%
% This function applies a cosinusoidal amplitude modulation of given rate
% and depth to an input sound. It returns the amplitude modulated sound.
% The function works either with monophonic and stereophonic sounds.
%
% SF: sample frequency of the tone in Hz
% SOUNDIN: the vector with the input sound
% MRATE: the modultation rate (in Hz)
% MDEPTH: the modulation depth. no modulation = 0; max modulation = 1
% MPHASE: starting phase of the modulator in radians (default is 0)
% ITD: ITD of the modulator in milliseconds. For positive values the signal
% comes from the left.


numberofsamples = length(soundin);
t = ((0:numberofsamples-1) / sf)';

% apply the modulation
if size(soundin, 2)==2
    
    modulatorL = zeros(numberofsamples, 1);
    modulatorR = modulatorL;
    
    ITD = ITD/10^3; % changes ITD value from milliseconds to seconds
    IPD = 2*pi*mrate*ITD;
    phaseL = mphase + IPD;
    phaseR = mphase;

    % generate and add each frequency component
    modulatorL = (1+cos(2*pi*mrate*t+phaseL))/2;
    modulatorL = 1-(modulatorL * mdepth);

    % generate and add each frequency component
    modulatorR = (1+cos(2*pi*mrate*t+phaseR))/2;
    modulatorR = 1-(modulatorR * mdepth);

    soundout = [soundin(:, 1).*modulatorL, soundin(:, 2).*modulatorR];
    % amplitude normalisation
    soundout = soundout / max(max(abs(soundout)));
else
    if nargin<5, mphase = 0; end

    modulator = zeros(numberofsamples, 1);
    % generate and add each frequency component
    modulator = (1+cos(2*pi*mrate*t+mphase))/2;
    modulator = 1-(modulator * mdepth);
    % apply the modulation
    soundout = soundin .* modulator;
    % amplitude normalisation
    soundout = soundout / max(abs(soundout));
end
% little correction
soundout = soundout * .999;