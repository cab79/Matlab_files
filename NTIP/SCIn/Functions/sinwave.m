function h=sinwave(h)%(A,freq,phase,samplerate,active,figureon)

% setup to: 
%   output either a single trials (h.Settings.trialdur>0) or multiple trials concatenated (h.Settings.trialdur=0)
%   vary across trials by intensity, duration, pitch or channel ('signal').
%   - each trial is create separately
%   vary within each trial by a pattern defined by intensity, duration, pitch or channel ('pattern')
%   - each pattern is created by a masking/splicing procedure

% https://uk.mathworks.com/help/matlab/import_export/record-and-play-audio.html#bsdl2eo-1
% "sound" simply plays the sound, "audioplayer" creates an object that
% allows puase and resume - maybe faster, but timing less important for
% fMRI.

% select channels
if isfield(h.Settings,'stimchan')
    chan = h.Settings.stimchan;
else
    chan = 1:h.Settings.nrchannels; % use all channels by default
end

% if ADAPTIVE
adaptive=0;
if isfield(h.Settings,'adaptive')
    if ~isempty(h.Settings.adaptive)
        adaptive=1;
    end
end

% if ODDBALL
oddball=0;
if isfield(h.Settings,'oddball')
    if ~isempty(h.Settings.oddball)
        oddball=1;
    end
end

% find trial(s) for which to create sin wave
if isfield(h,'i') % if a single trial has been defined in "trials" design
    trials = h.i;
else
    trials = 1:length(h.Seq.signal); % otherwise concatenate all trials
end

% trial loop
h.Seq.stimseq = [];
% sample ending each trial
h.Seq.trialend = [];
% instantaneous phase at end of trial
if ~isfield(h,'iphase')
    h.iphase = [];
end
for tr = trials
    
    % create temporary variables for intensity, pitch and total duration
    h.inten = h.Settings.inten;
    h.pitch = h.Settings.f0;
    % if calculating all trials, require totdur:
    if h.Settings.trialdur==0 && ~isfield(h,'totdur') % Use calculated duration by default.
        if isfield(h.Settings,'totdur')
            h.totdur = h.Settings.totdur; 
        end
    end
    % duration of each stimulus
    h.dur = h.Settings.stimdur; 

    % oddballs: select intensity and pitch
    if oddball
        if ~adaptive
            if strcmp(h.Settings.oddball,'channel')
                chan = chan(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddball,'intensity')
                h.inten = h.inten(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddball,'pitch')
                h.pitch = h.pitch(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddball,'duration')
                h.dur = h.dur(h.Seq.signal(tr),:);
            end
        elseif adaptive
            if isfield(h,'s')
                varlevel = h.s.StimulusLevel;
            else
                varlevel = h.Settings.adaptive.startinglevel;
            end
            if strcmp(h.Settings.oddball,'pitch')
                h.pitch = [h.Settings.f0(1), (h.Settings.f0(1)+varlevel)]; % create new pitch pair
                h.pitch = h.pitch(h.Seq.signal(tr));
            end
        end
    end
    
    % time index of sound
    t = transpose((1:sum(h.dur)*h.Settings.fs)/h.Settings.fs);
    
    % phase to add to stim
    try
        phadd = h.iphase(:,tr-1); % iphase at end of previous stim
    catch
        phadd = zeros(length(chan),1);
    end
        
    % construct waveform and its alternates to create patterns by masking/splicing
    h.mwav=zeros(h.Settings.nrchannels,length(t));
    if isfield(h.Settings,'df') && length(chan)==2 % if pitch is different in two channels
        df=1;
    else
        df=0;
    end
       
    if isfield(h.Settings,'pattern')
        if strcmp(h.Settings.pattern,'pitch') % pitch changes
            if length(h.dur)~=length(h.pitch)
                error('num column of stimdur must equate to number of pitches');
            end
            for i = 1:length(h.pitch)
                h.mwav(chan(1),:,i) = h.inten(1) *sin(2*pi*(h.pitch(i))*t + phadd(1) + 2*pi*h.pitch(i)/h.Settings.fs);
                if df; h.mwav(chan(2),:,i) = h.inten(1) *sin(2*pi*(h.pitch(i)+h.Settings.df)*t + phadd(2) + 2*pi*(h.pitch(i)+h.Settings.df)/h.Settings.fs);end
            end
        elseif strcmp(h.Settings.pattern,'intensity') % intensity changes
            if length(h.dur)~=length(h.inten)
                error('num column of stimdur must equate to number of intensities');
            end
            for i = 1:length(h.inten)
                h.mwav(chan(1),:,i) = h.inten(i) *sin(2*pi*(h.pitch(1))*t + phadd(1) + 2*pi*h.pitch(1)/h.Settings.fs);
                if df; h.mwav(chan(2),:,i) = h.inten(i) *sin(2*pi*(h.pitch(1)+h.Settings.df)*t + phadd(2) + 2*pi*(h.pitch(1)+h.Settings.df)/h.Settings.fs);end
            end
        end
    else
        % construct the player object: left
        h.mwav(chan(1),:) = h.inten(1) *sin(2*pi*h.pitch(1)*t + phadd(1) + 2*pi*h.pitch(1)/h.Settings.fs); % plus phaseshift plus increment
        % construct the player object: right
        if df; h.mwav(chan(2),:) = h.inten(1) *sin(2*pi*(h.pitch(1)+h.Settings.df)*t + phadd(2) + 2*pi*(h.pitch(1)+h.Settings.df)/h.Settings.fs);end
    end
    
    % if only one channel has so far been defined:
    if df==0
        h.mwav(chan,:,:) = repmat(h.mwav(chan(1),:,:),length(chan),1); 
    end
    
    %else
    %    wav = h.inten(1) *sin(2*pi*h.pitch(1)*t + phadd + 2*pi*h.pitch(1)/h.Settings.fs);

        % alternate sin waves
        %if h.Settings.fpitch>0 % pitch changes
        %    wav(:,2) = h.inten(1) *sin(2*pi*(h.pitch(2))*t);
        %elseif h.Settings.finten>0 % intensity changes
        %    wav(:,2) = h.inten(2) *sin(2*pi*(h.pitch(1))*t);
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
    if ~adaptive 
        h.inten_atten = h.Settings.atten; 
    elseif adaptive && oddball
        if strcmp(h.Settings.oddball,'intensity')
            h.inten_atten = [h.Settings.atten, (h.Settings.atten+varlevel)];
            h.inten_atten = h.inten_atten(h.Seq.signal(tr));
        else
            h.inten_atten = h.Settings.atten; 
        end
    end
    % find max rms_sound_dB
    h.mwav_orig = h.mwav; % non-attenuated version
    temp = reshape(permute(h.mwav,[2,1,3]),length(t),[]);
    for i = 1:size(temp,2)
        rms_sound_dB(i) = norm(temp(:,i))/sqrt(length(temp(:,i)));
    end
    rms_sound_dB = max(rms_sound_dB);
    ratio = min(1,(10^(h.inten_atten/20))/rms_sound_dB); % should always be smaller than 1
    h.mwav = ratio*h.mwav;
    
    if size(h.mwav,3)>1
        if isfield(h.Settings,'stimdurtype')
            if strcmp(h.Settings.stimdurtype,'rand')
                h.dur = h.dur(randperm(length(h.dur)));
            end
        end
        % create mask for temporal patterns
        mask = [];
        for i = 1:size(h.mwav,3)
            mask = [mask; i*ones(round(h.dur(i)*h.Settings.fs),1)]; 
        end
        % splice them together
        tempall = [];
        tempall2 = [];
        for i = 1:size(h.mwav,3)
            temp = h.mwav(:,mask==i,i);
            temp2 = h.mwav_orig(:,mask==i,i);
            % apply tapering
            if isfield(h.Settings,'Tukey')
                if h.Settings.Tukey>0
                    temp = temp.*repmat(tukeywin(size(temp,2),h.Settings.Tukey)',size(temp,1),1);
                end
            end
            tempall(:,mask==i) = temp;
            tempall2(:,mask==i) = temp2;
        end
        h.mwav=tempall;
        h.mwav_orig=tempall2;
        
    end
        
    % instantaneous phase at end of stim
    h.iphase(:,tr) = asin(h.mwav_orig(:,end));
    
    if length(trials)>1
        disp(['CREATE SEQUENCE: Trial ' num2str(tr) '/' num2str(length(trials)) ' appended']);
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
    h.Seq.stimseq = [h.Seq.stimseq h.mwav];
    h.Seq.trialend(tr,1) = size(h.Seq.stimseq,2);
end

function h = sinwave_cont(h)

t = transpose((1:h.totdur*h.Settings.fs)/h.Settings.fs);

% construct the player object: left
x = sin(2*pi*h.Settings.f0*t);
% construct the player object: right
y = sin(2*pi*(h.Settings.f0+h.Settings.df)*t);

if length(i12)~=length(t)
    error('lengths do not match')
end

% define durations - for duropt durations
i12 = [];
to = [];
for i = 1:length(randind)
    i12 = [
        i12;
        ones(round(h.Settings.duropt(randind(i),1)*h.Settings.fs),1); 
        zeros(round(h.Settings.duropt(randind(i),2)*h.Settings.fs),1)];
    if h.Settings.tone_dur>0
        to = [
            to;
            ones(round(h.Settings.tone_dur*h.Settings.fs),1); zeros(round((h.Settings.duropt(randind(i),1)-h.Settings.tone_dur)*h.Settings.fs),1);
            ones(round(h.Settings.tone_dur*h.Settings.fs),1); zeros(round((h.Settings.duropt(randind(i),2)-h.Settings.tone_dur)*h.Settings.fs),1);
        ];
    else
        to = ones(length(t),1);
    end
    disp(['SETUP SEQUENCE: Pair ' num2str(i) '/' num2str(length(randind)) ' appended']);
end

if length(i12)~=length(t)
    error('lengths do not match')
end

% pitch changes
if h.Settings.fpitch>0
    % alternate sin waves
    x2 = sin(2*pi*(h.Settings.f0+h.Settings.pitchdiff)*t);
    y2 = sin(2*pi*(h.Settings.f0+h.Settings.df+h.Settings.pitchdiff)*t);

    % define durations - for identical durations
    %i1 = ones((1/fpitch)*h.Settings.fs,1); 
    %i2 = zeros((1/fpitch)*h.Settings.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(h.Settings.fs/(length(i1)+length(i2))),1);

    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end

% intensity changes
if h.Settings.finten>0
    % alternate sin waves
    x2 = h.Settings.intendiff*sin(2*pi*(h.Settings.f0)*t);
    y2 = h.Settings.intendiff*sin(2*pi*(h.Settings.f0+h.Settings.df)*t);

    % define durations - for identical durations
    %i1 = ones((1/finten)*h.Settings.fs,1); 
    %i2 = zeros((1/finten)*h.Settings.fs,1); 
    %i12 = repmat([i1;i2],mod_dur*(h.Settings.fs/(length(i1)+length(i2))),1);

    % splice them in
    x(find(i12)) = x2(find(i12));
    y(find(i12)) = y2(find(i12));
end
x = x.*to;
y = y.*to;

h.Seq.signal = audioplayer([x y], h.Settings.fs);
h.Seq.blocks = ones(1,1:length(x));