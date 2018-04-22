function h = contruct_wave(h)
% CONSTRUCT WAVEFORM

% set default to sin wave
if ~isfield(h.Settings,'wavetype')
    h.Settings.wavetype = 'sin';
elseif isempty(h.Settings.wavetype)
    h.Settings.wavetype = 'sin';
end

if ~isfield(h,'chan')
    h.chan = 1;
end
if ~isfield(h,'trialtype')
    h.trialtype.freqpattern = 0;
    h.trialtype.intenpattern = 0;
end
if ~isfield(h,'freq')
    h.freq = h.Settings.stim(h.trialstimnum).f0;
end

if isfield(h.Settings,'df') && length(h.chan)==2 % if freq is different in two channels
    df=1;
    if isfield(h,'entrainfreq')
        try
            df_freq = str2double(get(h.entrainfreq,'string'));
        catch
            df_freq = str2double(h.entrainfreq);
        end
        if df_freq==0 || isnan(df_freq)
            df_freq=h.Settings.stim(h.trialstimnum).df;
        end
    else
        df_freq=h.Settings.stim(h.trialstimnum).df;
    end
else
    df=0;
end

% randomise order?
rs = 1:length(h.dur);
if ~isfield(h,'resp_probe')
    h.resp_probe=0;
end
if isfield(h.Settings.stim(h.trialstimnum),'durtype') && ~h.resp_probe
    if strcmp(h.Settings.stim(h.trialstimnum).durtype,'rand')
        rs = randperm(length(h.Settings.stimrandind));
        h.dur(h.Settings.stimrandind) = h.dur(h.Settings.stimrandind(rs));
        if strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'intensity')
            h.inten(h.Settings.stimrandind) = h.inten(h.Settings.stimrandind(rs));
        elseif strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'pitch') || strcmp(h.Settings.stim(h.trialstimnum).patternmethod,'freq')
            h.freq(h.Settings.stimrandind) = h.freq(h.Settings.stimrandind(rs));
        end
    end
end

h.mwav=[];
% for each stimdur
for i = 1:length(h.dur)

    % time index of sound
    t{i} = transpose((1:h.dur(i)*h.Settings.fs)/h.Settings.fs);

    % start phase to add to stim
    if h.alignphase
        try
            if i==1 % use previous trial's end phase
                phadd = h.iphase(:,tr-1,1,end); % phase at end of previous stim
                phadd = phadd+(h.iphase(:,tr-1,2,end)<0).*(pi-2*h.iphase(:,tr-1,1,end)); % adjust when slope of sinwave is negative, which gives wrong phase values
            else
                phadd = h.iphase(:,tr,1,i-1); % phase at end of previous stim
                phadd = phadd+(h.iphase(:,tr,2,i-1)<0).*(pi-2*h.iphase(:,tr,1,i-1)); % adjust when slope of sinwave is negative, which gives wrong phase values    
            end
        catch
            phadd = zeros(length(h.chan),1);
        end
    else
        phadd = zeros(length(h.chan),1);
    end

    % initialise
    mwav{i}=zeros(h.Settings.stim(h.trialstimnum).nrchannels,length(t{i}));
    temp_sin{i}=zeros(h.Settings.stim(h.trialstimnum).nrchannels,length(t{i}));

    % pitch/inten-specific 
    if h.trialtype.freqpattern
        %h.freq = h.freq(rs); % randomise?
        if length(h.dur)~=length(h.freq)
            error('num column of stimdur must equate to number of frequencies');
        end
        usefreq = h.freq(i);
        useinten = h.inten(1);
    elseif h.trialtype.intenpattern
        %h.inten = h.inten(rs); % randomise?
        if length(h.dur)~=length(h.inten)
            error('num column of stimdur must equate to number of intensities');
        end
        usefreq = h.freq(1);
        useinten = h.inten(i);
    else % no pattern
        usefreq = h.freq(i);
        useinten = h.inten(i);
    end

    % construction of mwav{i}
    if strcmp(h.Settings.wavetype,'sin')
        
        % are we using decibels? If so, set to 1 so that volume is adjusted
        % later
        if strcmp(h.Settings.stim(h.trialstimnum).inten_type,'dB') % decibels scale
            useinten=1;
        end
        
        mwav{i}(h.chan(1),:) = useinten *sin(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs);
        if df; mwav{i}(h.chan(2),:) = useinten *sin(2*pi*(usefreq+df_freq)*t{i} + phadd(2) + 2*pi*(usefreq+df_freq)/h.Settings.fs);end
        temp_sin{i} = mwav{i};
    elseif strcmp(h.Settings.wavetype,'square')
        if isfield(h,'dutycycle')
            if h.dutycycle>0
                useinten = str2double(h.dutycycle);
            end
        end
        mwav{i}(h.chan(1),:) = square(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs,useinten);
        temp_sin{i}(h.chan(1),:) = sin(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs); % for phase estimation
    elseif strcmp(h.Settings.wavetype,'step')
        %h.Settings.fs = 100;
        %usefreq = 1; % rise-fall rate /s
        %useinten = 3; % volts
        %h.dur = 3; % duration at max intensity
        %h.Settings.trialdur = 10;

        rise = transpose((1:(useinten/usefreq)*h.Settings.fs)/h.Settings.fs);
        fall = flip(rise);
        peak = useinten*ones(h.dur*h.Settings.fs,1);
        trough = zeros((h.Settings.trialdur*h.Settings.fs)-length(rise)-length(fall)-length(peak),1);

        mwav{i} = [rise;peak;fall;trough]';
    end

    % instantaneous phase and direction at end of stim
    h.iphase(:,h.tr,:,i) = [asin(temp_sin{i}(:,end-1)), temp_sin{i}(:,end)-temp_sin{i}(:,end-1)];% mwav_orig(:,end-1)-mwav_orig(:,end)];

    %optional plots
    %disp_phase=h.iphase(:,tr,:,i)
    %if i>1
    %    close all
    %    figure;
    %    subplot(1,3,1);
    %    plot([mwav{i-1}(1,end-500:end),mwav{i}(1,1:500)]);
    %    subplot(1,3,2);
    %    plot([mwav{i-1}(2,end-500:end),mwav{i}(2,1:500)]);
    %    subplot(1,3,3);
    %    plot([mwav{i-1}(2,end-500:end)-mwav{i-1}(1,end-500:end),mwav{i}(2,1:500)-mwav{i}(1,1:500)]);
    %end

    % apply tapering on each part of the stim pattern?
    if isfield(h.Settings.stim(h.trialstimnum),'Tukey')
        if h.Settings.stim(h.trialstimnum).Tukeytype==1
            mwav{i} = mwav{i}.*repmat(tukeywin(size(mwav{i},2),h.Settings.stim(h.trialstimnum).Tukey)',size(mwav{i},1),1);
        end
    end

    % concatenate
    h.mwav = [h.mwav mwav{i}];

end

% create tsum
tsum = sum(cellfun(@length,t));


% construct waveform and its alternates to create patterns by masking/splicing
%h.mwav=zeros(h.Settings.nrchannels,length(t));
%if isfield(h.Settings,'df') && length(chan)==2 % if pitch is different in two channels
%    df=1;
%else
%    df=0;
%end

% if freqpattern
%     if length(h.dur)~=length(h.freq)
%         error('num column of stimdur must equate to number of pitches');
%     end
%     for i = 1:length(h.freq)
%         h.mwav(chan(1),:,i) = h.inten(1) *sin(2*pi*(h.freq(i))*t + phadd(1) + 2*pi*h.freq(i)/h.Settings.fs);
%         if df; h.mwav(chan(2),:,i) = h.inten(1) *sin(2*pi*(h.freq(i)+df_freq)*t + phadd(2) + 2*pi*(h.freq(i)+df_freq)/h.Settings.fs);end
%     end
% elseif intenpattern
%     if length(h.dur)~=length(h.inten)
%         error('num column of stimdur must equate to number of intensities');
%     end
%     for i = 1:length(h.inten)
%         h.mwav(chan(1),:,i) = h.inten(i) *sin(2*pi*(h.freq(1))*t + phadd(1) + 2*pi*h.freq(1)/h.Settings.fs);
%         if df; h.mwav(chan(2),:,i) = h.inten(i) *sin(2*pi*(h.freq(1)+df_freq)*t + phadd(2) + 2*pi*(h.freq(1)+df_freq)/h.Settings.fs);end
%     end
% end


% otherwise, constuct waveform
% if ~any(h.mwav(:))
%     % construct the player object: left
%     h.mwav(chan(1),:) = h.inten(1) *sin(2*pi*h.freq(1)*t + phadd(1) + 2*pi*h.freq(1)/h.Settings.fs); % plus phaseshift plus increment
%     % construct the player object: right
%     if df; h.mwav(chan(2),:) = h.inten(1) *sin(2*pi*(h.freq(1)+df_freq)*t + phadd(2) + 2*pi*(h.freq(1)+df_freq)/h.Settings.fs);end
% end
mono=0;
if isfield(h.Settings,'monostereo')
    if strcmp(h.Settings.monostereo,'mono')
        mono=1;
    end
end

% if only one channel has so far been defined:
if df==0 || mono
    h.mwav(h.chan,:,:) = repmat(h.mwav(h.chan(1),:,:),length(h.chan),1); 
end
% monaural beats
if isfield(h.Settings,'monaural')
    if h.Settings.monaural && size(h.mwav,1)==2
        h.mon=h.mwav;
        h.mon(h.chan(1),:,:) = h.mon(h.chan(1),:,:) - h.mon(h.chan(2),:,:);
        h.mon(h.chan(2),:,:) = h.mon(h.chan(1),:,:);
    end
end

%else
%    wav = h.inten(1) *sin(2*pi*h.freq(1)*t + phadd + 2*pi*h.freq(1)/h.Settings.fs);

    % alternate sin waves
    %if h.Settings.fpitch>0 % pitch changes
    %    wav(:,2) = h.inten(1) *sin(2*pi*(h.freq(2))*t);
    %elseif h.Settings.finten>0 % intensity changes
    %    wav(:,2) = h.inten(2) *sin(2*pi*(h.freq(1))*t);
    %end

%    h.mwav(chan,:) = repmat(wav,length(chan),1); 
%end


% amplitude normalisation
%wav = wav' / max(max(abs(wav)));
% little correction
%wav = wav * .999;

% This attenuates the overall level of a sound of a certain number
% of decibels. The reference for the attenuation is 0-dB FS (i.e., full
% scale). The function attenuates the sound if the input number of decibels
% is negative and amplify the sound otherwise.
% If the sound is stereophonic the desired attenuation (att_dB) is imposed
% on the loudest channel.
% Please, note that a sinusoidal signal with peak amplitude at 1 has an
% overall level of -3.02 dB.
% ATT_DB: attenuation/amplification in dB
% DIGITAL SOUNDS HAVE A MAXIMUM DYNAMIC RANGE OF 96 dB FOR AT 16 BITS
% RESOLUTION:
%       20*log10(2^16)=96.33
if ~isfield(h,'varlevel')
    h.varlevel=0;
end
inten_atten=0;
if isfield(h.Settings.stim(h.trialstimnum),'attenchan')
    if any(ismember(h.chan,h.Settings.stim(h.trialstimnum).attenchan))
        
        % overall attenuation from GUI
        if isfield(h,'vol_atten')
            try
                inten_atten = str2double(get(h.vol_atten,'string'));
            catch
                inten_atten = str2double(h.vol_atten);
            end
        end
        
        % overall attentuation from settings/sequence
        if iscell(h.Settings.stim(h.trialstimnum).atten)
            if strcmp(h.Settings.stim(h.trialstimnum).atten{1},'inten_diff')
                if ~h.seqtype.adapt && ~h.seqtype.thresh
                    inten_atten = inten_atten + h.inten_diff * h.Settings.stim(h.trialstimnum).atten{2};
                else
                    if isfield(h,'s')
                        inten_atten = h.s.a(strcmp({h.Settings.adaptive(:).type}, 'discrim')).StimulusLevel;
                    else
                        if h.seqtype.adapt
                            inten_atten = h.Settings.adaptive.startinglevel;
                        else
                            inten_atten = h.Settings.threshold.startinglevel;
                        end  
                    end
                end
            end
        else % use numeric value from settings
            inten_atten = inten_atten + h.Settings.stim(h.trialstimnum).atten; 
        end
        
        % stim-specific attentuation from sequence
        if ~h.seqtype.adapt && ~h.seqtype.thresh
            h.inten_atten = inten_atten+h.varlevel;
        elseif (h.seqtype.adapt || h.seqtype.thresh) && h.seqtype.oddball && strcmp(h.Settings.stim(h.trialstimnum).inten_type,'dB')
            if strcmp(h.Settings.oddballmethod,'intensity')
                h.inten_atten = [inten_atten-h.varlevel/2, (inten_atten+h.varlevel/2)];
                h.inten_atten = h.inten_atten(h.Seq.signal(h.trialstimnum,h.tr));
            else
                h.inten_atten = inten_atten+h.varlevel;
            end
        elseif h.seqtype.thresh && strcmp(h.Settings.stim(h.trialstimnum).inten_type,'dB')
            if strcmp(h.Settings.threshold.type,'intensity')
                h.inten_atten = inten_atten+h.varlevel;
            end
        end
    else
        h.inten_atten = 0; 
    end
else
    h.inten_atten = 0; 
end

if h.inten_atten
    h.mwav = attenute_sound(h.mwav,h.inten_atten);
end

if strcmp(h.Settings.stim(h.trialstimnum).inten_type,'dB')
    h.inten_out = h.inten+h.inten_atten;
else
    h.inten_out = h.inten;
end


if isfield(h.Settings.stim(h.trialstimnum),'Tukey')
    if h.Settings.stim(h.trialstimnum).Tukeytype==2
        h.mwav = h.mwav.*repmat(tukeywin(size(h.mwav,2),h.Settings.stim(h.trialstimnum).Tukey)',size(h.mwav,1),1);
        if isfield(h,'mon');h.mon = h.mon.*repmat(tukeywin(size(h.mon,2),h.Settings.stim(h.trialstimnum).Tukey)',size(h.mon,1),1);end
    end
end


if length(h.trials)>1
    disp(['CREATE SEQUENCE: Trial ' num2str(h.tr) '/' num2str(length(h.trials)) ' appended']);
end

% THE FOLLOWING PATTERNS WILL BE DEFINED WITHIN SETTINGS
% generate times within active period
%xmax=h.dur;
%n=h.Settings.npulses;

%if n>1
%    if strcmp(h.Settings.sinwave,'rand');
%        x=sort(rand(1,n*2)*xmax);
%    elseif strcmp(h.Settings.sinwave,'reg');
%        x=sort([1/h.Settings.fs 1/(n*2):1/(n*2):(1-1/(n*2))]*xmax);
%    end
%else
%    x=sort([1/h.Settings.fs,1]*xmax);
%end

%x=reshape(x,2,length(x)/2);

%add waveform into those random times
%sig = [repmat(wav(1:end-1),1,n) wav(end)];
%h.wav = zeros(size(sig));

%for i = 1:size(x,2)
%    temp = sig(max(1,round(x(1,i)*h.Settings.fs)):round(x(2,i)*h.Settings.fs));
%    if isfield(h.Settings,'Tukey')
%        if h.Settings.Tukey>0
%            temp = temp.*tukeywin(length(temp),h.Settings.Tukey)';
%        end
%    end
%    h.wav(max(1,round(x(1,i)*h.Settings.fs)):round(x(2,i)*h.Settings.fs)) = temp;

%    ha = area(x(:,i)', [1 1]);
%    axis([0 active 0 1])
%    hold on
%end


% display figure
%if figureon
%    figure
%    plot(0:1/samplerate:active,sigB)
%end
if isfield(h,'mon')
    h.Seq.stimseq = [h.Seq.stimseq h.mon];
else
    h.Seq.stimseq = [h.Seq.stimseq h.mwav];
end

if strcmp(h.Settings.design,'continuous')
    if isfield(h,'i') % add to trialend according to true position in sequence
        if h.i>1
            h.Seq.trialend(h.i+h.Settings.ntrialsahead-1,1) = h.Seq.trialend(h.i-1+h.Settings.ntrialsahead-1,1)+size(h.mwav,2);
            h.Seq.PatternSamples{h.i+h.Settings.ntrialsahead-1,1} = cumsum(h.dur*h.Settings.fs);
        elseif h.i==1
            h.Seq.trialend(h.i,1) = size(h.mwav,2);
            h.Seq.PatternSamples{h.i,1} = cumsum(h.dur*h.Settings.fs);
        end
    else % experiment not started, so add to trialend using tr (which will start from 1)
        try
            h.Seq.trialend(h.tr,1) = h.Seq.trialend(h.tr-1,1)+size(h.mwav,2);
        catch
            h.Seq.trialend(h.tr,1) = size(h.mwav,2);
        end
        h.Seq.PatternSamples{h.tr,1} = cumsum(h.dur*h.Settings.fs);
    end 
end

%figure;plot(h.Seq.stimseq(2,h.Seq.PatternSamples{tr,1}(1)-500:h.Seq.PatternSamples{tr,1}(1)+500));