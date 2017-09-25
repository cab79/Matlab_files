function noise = GenerateNoise (sf, dur, noisetype, lf, hf)
%
% noise = GenerateNoise (sf, dur)
% filterednoise = GenerateNoise (sf, dur, noisetype, f1)
% filterednoise = GenerateNoise (sf, dur, noisetype, f1, f2)
% 
% This function generates a Gaussian noise either white (default option if
% noisetype is omitted), lowpass, highpass, bandpass or notched (i.e.,
% bandstop).
%
% SF: sample frequency of the noise in Hz
% DUR: noise duration in msec
% NOISETYPE: a string with the kind of noise has to be generated. 
% Accepted noisetype are 'lowpass', 'highpass', 'bandpass' and 'notched'
% F1: lowest frequency of the filter in Hz
% F2: highest frequency of the filter in Hz. This value can be omitted for
% highpass and lowpass filtered noises.
%
% % EXAMPLE: generate a lowpas filtered noise of 500-ms with cutoff
% frequency of 5000-Hz
% filterednoise = GenerateFilteredNoise (44100, 500, 'lowpass', 5000);

if nargin < 3
    % set general variables
    nf = sf/2;                              % Nyquist frequency
    numberofsamples = round(sf*(dur/1000)); % number of samples
    nh = round(numberofsamples/2);                 % half number of samples

    % make noise
    rand('state',sum(100 * clock));     % initialize random seed
    noise = randn(numberofsamples, 1);  % Gaussian noise
    noise = noise / max(abs(noise));    % -1 to 1 normalization
    % amplitude normalization
    noise = noise/max(abs(noise));
    noise = noise*.999;
else
    if nargin<5, hf=lf; end;
    if hf<lf
        error('hf cannot be greater than than lf');
    end;

    % set general variables
    nf = sf/2;                                  % Nyquist frequency
    numberofsamples = round(sf*(dur/1000));     % number of samples
    nh = round(numberofsamples/2);              % half number of samples
    % make noise
    rand('state',sum(100 * clock));             % initialize random seed
    noise = randn(numberofsamples, 1);          % Gaussian noise
    noise = noise / max(abs(noise));            % -1 to 1 normalization
    % =========================================================================
    % set variables for filter
    lp = lf * (dur/1000); % lf point in frequency domain
    hp = hf * (dur/1000); % hf point in frequency domain
    % design filter
    switch noisetype
        case 'lowpass'
            filter = zeros(numberofsamples, 1);
            filter(1 : lp, 1) = 1;
            filter(numberofsamples - lp : numberofsamples, 1) = 1;
        case 'highpass'
            filter = ones(numberofsamples, 1);
            filter(1 : hp, 1) = 0;
            filter(numberofsamples - hp : numberofsamples, 1) = 0;
        case 'bandpass'
            filter = zeros(numberofsamples, 1);
            filter(lp : hp, 1) = 1;
            filter(numberofsamples - hp : numberofsamples - lp, 1) = 1;
        case 'notched'
            filter = ones(numberofsamples, 1);
            filter(lp : hp, 1) = 0;
            filter(numberofsamples - hp : numberofsamples - lp, 1) = 0;
        otherwise
            error('unknown noisetype');
    end;
    % =========================================================================
    % do filter
    noise = fft(noise);                 % FFT
    noise = noise .* filter;            % filtering
    noise = ifft(noise);                % inverse FFT
    noise = real(noise);
    % amplitude normalization
    noise = noise/max(abs(noise));
    noise = noise*.999;
end;