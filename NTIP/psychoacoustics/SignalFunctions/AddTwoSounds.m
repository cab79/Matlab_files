function soundout = AddTwoSounds(sf, soundin1, soundin2, onset)
%
% soundout = AddTwoSounds(sf, soundin1, soundin2, soundin2onset)
%
% This function adds pairs of sounds. It receives two input sounds and
% returns a single sound. The second sound will be switched on at onset
% time, relatively to the beginning of the first sound. The function works
% also if sounds do not to overlap in time, for example, when the second
% sound is switched well after the end of the first sound. If you need to
% strictly concatenate the sounds (i.e., no silent gap between first and
% second sound) use the ConcatenateSounds function instead. AddTwoSounds
% works with either monophonic and stereophonic sounds.
% 
% ATTENTION: the function does not normalize the amplitude of the resulting
% sound. Therefore, the resulting sound can be clipped.
%
% SF: sample frequency in Hz
% SOUNDIN1: the sound vector of the first sound
% SOUNDIN2: the sound vector of the second sound
% SOUNDIN2ONSET: onset of the second sound relatively to the first (in ms)
%
% % EXAMPLE 1: create the stimulus for a simultaneous masking experiment.
% % A brief sine tone will be added to a band of bandpass noise. 
% sf = 44100;
% noise = GenerateNoise(sf, 300, 'bandpass', 600, 1400);
% noise = GenerateEnvelope(sf, noise);
% noise = AttenuateSound(noise, -40);
% tone = GenerateTone(sf, 1000, 20);
% tone = GenerateEnvelope(sf, tone);
% tone = AttenuateSound(tone, -10);
% TonePlusNoise = AddTwoSounds(sf, noise, tone, 140);
% sound(TonePlusNoise, sf)
  
if size(soundin1, 2)==2
    silencepre = zeros(round((onset/1000)*sf), 2);
else
    silencepre = zeros(round((onset/1000)*sf), 1);
end

if (length(silencepre)+length(soundin2)) < length(soundin1) && size(soundin1, 2)==2
    silencepost = zeros(length(soundin1)-(length(silencepre)+length(soundin2)), 2);
    soundin2 = [silencepre; soundin2; silencepost];
    soundout = soundin1 + soundin2;
elseif (length(silencepre)+length(soundin2)) > length(soundin1) && size(soundin1, 2)==2
    silencepost = zeros((length(silencepre)+length(soundin2))-length(soundin1), 2);
    soundin1 = [soundin1; silencepost];
    soundin2 = [silencepre; soundin2];
    soundout = soundin1 + soundin2;
elseif (length(silencepre)+length(soundin2)) == length(soundin1) && size(soundin1, 2)==2
    soundin2 = [silencepre; soundin2];
    soundout = soundin1 + soundin2;
elseif (length(silencepre)+length(soundin2)) < length(soundin1) && size(soundin1, 2)==1
    silencepost = zeros(length(soundin1)-(length(silencepre)+length(soundin2)), 1);
    soundin2 = [silencepre; soundin2; silencepost];
    soundout = soundin1 + soundin2;
elseif (length(silencepre)+length(soundin2)) > length(soundin1) && size(soundin1, 2)==1
    silencepost = zeros((length(silencepre)+length(soundin2))-length(soundin1), 1);
    soundin1 = [soundin1; silencepost];
    soundin2 = [silencepre; soundin2];
    soundout = soundin1 + soundin2;
elseif (length(silencepre)+length(soundin2)) == length(soundin1) && size(soundin1, 2)==1
    soundin2 = [silencepre; soundin2];
    soundout = soundin1 + soundin2;
end