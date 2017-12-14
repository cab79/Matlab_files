function y = attenute_sound(y,atten)

temp = reshape(permute(y,[2,1,3]),size(y,2),[]);
% find max rms_sound_dB
%h.mwav_orig = h.mwav; % non-attenuated version
for i = 1:size(temp,2)
    rms_sound_dB(i) = norm(temp(:,i))/sqrt(length(temp(:,i)));
end
rms_sound_dB = max(rms_sound_dB);
ratio = min(1,(10^(atten/20))/rms_sound_dB); % should always be smaller than 1
y = ratio*y;
end