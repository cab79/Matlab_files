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

dbstop if error

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
if isfield(h.Settings,'oddballmethod')
    if ~isempty(h.Settings.oddballmethod)
        oddball=1;
    end
end

% find trial(s) for which to create sin wave
if strcmp(h.Settings.design,'trials') && isfield(h,'i') % if a single trial has been defined in "trials" design
    trials = h.i;
elseif strcmp(h.Settings.design,'continuous')
    if h.Settings.ntrialsahead>0 && isfield(h,'i') % experiment is already running
        trials = h.i+(h.Settings.ntrialsahead-1);
    elseif h.Settings.ntrialsahead>0 && ~isfield(h,'i') % buffer needs pre-filling prior to experiment starting
        trials = 1:h.Settings.ntrialsahead;
    else
        trials = 1:length(h.Seq.signal); % otherwise concatenate all trials
    end
end 

% trial loop
h.Seq.stimseq = [];
% sample ending each trial
%h.Seq.trialend = [];
% instantaneous phase at end of trial
if ~isfield(h,'iphase')
    h.iphase = [];
end
for tr = trials
    % create temporary variables for intensity, pitch and duration
    h.inten = h.Settings.inten;
    h.pitch = h.Settings.f0;
    h.dur = h.Settings.stimdur; 
    % if calculating all trials, requires totdur:
    if strcmp(h.Settings.design,'continuous') && ~isfield(h,'totdur') % Use calculated duration by default.
        if isfield(h.Settings,'totdur')
            h.totdur = h.Settings.totdur; 
        end
    end

    % condition method
    if isfield(h.Settings,'conditionmethod')
        if ~isempty(h.Settings.conditionmethod)
            if iscell(h.Settings.conditionmethod)
                for i = 1:length(h.Settings.conditionmethod)
                    conditionmethod = h.Settings.conditionmethod{i};
                    if strcmp(conditionmethod,'pitch')
                        if iscell(h.Settings.conditionvalue)
                            h.pitch = h.Settings.conditionvalue{h.Seq.signal(tr),i};
                        else
                            h.pitch = h.Settings.conditionvalue(i,h.Seq.signal(tr));
                        end
                    end
                    if strcmp(conditionmethod,'intensity')
                        if iscell(h.Settings.conditionvalue)
                            h.inten = h.Settings.conditionvalue{h.Seq.signal(tr),i};
                        else
                            h.inten = h.Settings.conditionvalue(i,h.Seq.signal(tr));
                        end
                    end
                end
            else
                error('h.Settings.conditionmethod must be a cell');
            end
        end
    end
    
    % oddball method
    if oddball
        if ~adaptive
            if iscell(h.Settings.oddballvalue)
                if size(h.Settings.oddballvalue,1)==1
                    oddval = h.Settings.oddballvalue{h.Seq.signal(tr)};
                else
                    oddval = h.Settings.oddballvalue{h.Seq.signal(tr),:};
                end
            else
                oddval = h.Settings.oddballvalue(h.Seq.signal(tr),:);
            end
            if strcmp(h.Settings.oddballmethod,'channel')
                chan=oddval;
            elseif strcmp(h.Settings.oddballmethod,'intensity')
                h.inten = oddval;
            elseif strcmp(h.Settings.oddballmethod,'pitch')
                h.pitch = oddval;
            elseif strcmp(h.Settings.oddballmethod,'duration')
                h.dur = oddval;
            end
        elseif adaptive
            if isfield(h,'s')
                varlevel = h.s.StimulusLevel;
            else
                varlevel = h.Settings.adaptive.startinglevel;
            end
            if strcmp(h.Settings.oddballmethod,'pitch')
                h.pitch = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+varlevel)]; % create new pitch pair
                h.pitch = h.pitch(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddballmethod,'duration') && strcmp(h.Settings.patternmethod,'pitch')
                if iscell(h.Settings.oddballvalue)
                    h.dur = h.Settings.oddballvalue{h.Seq.signal(tr),:};
                else
                    h.dur = h.Settings.oddballvalue(h.Seq.signal(tr),:);
                end
                h.pitch = [h.Settings.patternvalue(1), (h.Settings.patternvalue(1)+varlevel)]; % create new pitch pair
            end
        end
    end
       
    %apply pitch pattern?
    pitchpattern=0;
    if isfield(h.Settings,'patternmethod')
        if strcmp(h.Settings.patternmethod,'pitch') % pitch changes
            pitchpattern=1;
            if ~(adaptive && strcmp(h.Settings.oddballmethod,'pitch')) && ~strcmp(conditionmethod,'pitch') % pitch already defined above in this case
                if isnumeric(h.Settings.patternvalue)
                    h.pitch = h.Settings.patternvalue;
                elseif iscell(h.Settings.patternvalue)
                    nDur = length(h.dur);
                    nPit = cellfun(@length,h.Settings.patternvalue);
                    h.pitch = h.Settings.patternvalue{nPit==nDur};
                end
            end
        end
    end
    % apply response probe?
    if isfield(h.Settings,'RPmethod')
        if strcmp(h.Settings.RPmethod,'pitch')
            if h.Seq.RP(tr)==1
                pitchpattern=1;
                h.pitch = h.Settings.RPvalue;
                h.dur = h.Settings.RPdur;
            end
        end
    end
    %apply intensity pattern?
    intenpattern=0;
    if isfield(h.Settings,'patternmethod')
        if strcmp(h.Settings.patternmethod,'intensity') % intensity changes
            intenpattern=1;
            if ~strcmp(conditionmethod,'intensity') % then already defined
                h.inten = h.Settings.patternvalue;
            end
        end
    end
    % apply response probe?
    resp_probe=0;
    if isfield(h.Settings,'RPmethod')
        if strcmp(h.Settings.RPmethod,'intensity')
            if h.Seq.RP(tr)==1
                resp_probe=1;
                h.inten = h.Settings.RPvalue;
                h.dur = h.Settings.RPdur;
            end
        end
    end
    
    %% CONSTRUCT WAVEFORM
    
    if isfield(h.Settings,'df') && length(chan)==2 % if pitch is different in two channels
        df=1;
        if isfield(h,'entrainfreq')
            df_freq = str2double(h.entrainfreq);
            if df_freq==0
                df_freq=h.Settings.df;
            end
        else
            df_freq=h.Settings.df;
        end
    else
        df=0;
    end
    
    % randomise order?
    rs = 1:length(h.dur);
    if isfield(h.Settings,'stimdurtype') && ~resp_probe
        if strcmp(h.Settings.stimdurtype,'rand')
            rs = randperm(length(h.Settings.stimrandind));
            h.dur(h.Settings.stimrandind) = h.dur(h.Settings.stimrandind(rs));
            if strcmp(h.Settings.patternmethod,'intensity')
                h.inten(h.Settings.stimrandind) = h.inten(h.Settings.stimrandind(rs));
            elseif strcmp(h.Settings.patternmethod,'pitch')
                h.pitch(h.Settings.stimrandind) = h.pitch(h.Settings.stimrandind(rs));
            end
        end
    end
    
    h.mwav=[];
    % for each stimdur
    for i = 1:length(h.dur)
        
        % time index of sound
        t{i} = transpose((1:h.dur(i)*h.Settings.fs)/h.Settings.fs);
        
        % start phase to add to stim
        try
            if i==1 % use previous trial's end phase
                phadd = h.iphase(:,tr-1,1,end); % phase at end of previous stim
                phadd = phadd+(h.iphase(:,tr-1,2,end)<0).*(pi-2*h.iphase(:,tr-1,1,end)); % adjust when slope of sinwave is negative, which gives wrong phase values
            else
                phadd = h.iphase(:,tr,1,i-1); % phase at end of previous stim
                phadd = phadd+(h.iphase(:,tr,2,i-1)<0).*(pi-2*h.iphase(:,tr,1,i-1)); % adjust when slope of sinwave is negative, which gives wrong phase values    
            end
        catch
            phadd = zeros(length(chan),1);
        end
        
        % initialise
        mwav{i}=zeros(h.Settings.nrchannels,length(t{i}));
        
        % pitch/inten-specific construction of mwav{i}
        if pitchpattern
            %h.pitch = h.pitch(rs); % randomise?
            if length(h.dur)~=length(h.pitch)
                error('num column of stimdur must equate to number of pitches');
            end
            mwav{i}(chan(1),:) = h.inten(1) *sin(2*pi*(h.pitch(i))*t{i} + phadd(1) + 2*pi*h.pitch(i)/h.Settings.fs);
            if df; mwav{i}(chan(2),:) = h.inten(1) *sin(2*pi*(h.pitch(i)+df_freq)*t{i} + phadd(2) + 2*pi*(h.pitch(i)+df_freq)/h.Settings.fs);end
        elseif intenpattern
            %h.inten = h.inten(rs); % randomise?
            if length(h.dur)~=length(h.inten)
                error('num column of stimdur must equate to number of intensities');
            end
            mwav{i}(chan(1),:) = h.inten(i) *sin(2*pi*(h.pitch(1))*t{i} + phadd(1) + 2*pi*h.pitch(1)/h.Settings.fs);
            if df; mwav{i}(chan(2),:) = h.inten(i) *sin(2*pi*(h.pitch(1)+df_freq)*t{i} + phadd(2) + 2*pi*(h.pitch(1)+df_freq)/h.Settings.fs);end
        else % no pattern
            mwav{i}(chan(1),:) = h.inten(i) *sin(2*pi*(h.pitch(i))*t{i} + phadd(1) + 2*pi*h.pitch(i)/h.Settings.fs);
            if df; mwav{i}(chan(2),:) = h.inten(i) *sin(2*pi*(h.pitch(i)+df_freq)*t{i} + phadd(2) + 2*pi*(h.pitch(i)+df_freq)/h.Settings.fs);end
        end
        
        % instantaneous phase and direction at end of stim
        h.iphase(:,tr,:,i) = [asin(mwav{i}(:,end-1)), mwav{i}(:,end)-mwav{i}(:,end-1)];% mwav_orig(:,end-1)-mwav_orig(:,end)];
        
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
        if isfield(h.Settings,'Tukey')
            if h.Settings.Tukeytype==1
                mwav{i} = mwav{i}.*repmat(tukeywin(size(mwav{i},2),h.Settings.Tukey)',size(mwav{i},1),1);
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
    
   % if pitchpattern
   %     if length(h.dur)~=length(h.pitch)
   %         error('num column of stimdur must equate to number of pitches');
   %     end
   %     for i = 1:length(h.pitch)
   %         h.mwav(chan(1),:,i) = h.inten(1) *sin(2*pi*(h.pitch(i))*t + phadd(1) + 2*pi*h.pitch(i)/h.Settings.fs);
   %         if df; h.mwav(chan(2),:,i) = h.inten(1) *sin(2*pi*(h.pitch(i)+df_freq)*t + phadd(2) + 2*pi*(h.pitch(i)+df_freq)/h.Settings.fs);end
   %     end
   % elseif intenpattern
   %     if length(h.dur)~=length(h.inten)
   %         error('num column of stimdur must equate to number of intensities');
   %     end
   %     for i = 1:length(h.inten)
   %         h.mwav(chan(1),:,i) = h.inten(i) *sin(2*pi*(h.pitch(1))*t + phadd(1) + 2*pi*h.pitch(1)/h.Settings.fs);
   %         if df; h.mwav(chan(2),:,i) = h.inten(i) *sin(2*pi*(h.pitch(1)+df_freq)*t + phadd(2) + 2*pi*(h.pitch(1)+df_freq)/h.Settings.fs);end
   %     end
   % end
    
    
    % otherwise, constuct waveform
   % if ~any(h.mwav(:))
   %     % construct the player object: left
   %     h.mwav(chan(1),:) = h.inten(1) *sin(2*pi*h.pitch(1)*t + phadd(1) + 2*pi*h.pitch(1)/h.Settings.fs); % plus phaseshift plus increment
   %     % construct the player object: right
   %     if df; h.mwav(chan(2),:) = h.inten(1) *sin(2*pi*(h.pitch(1)+df_freq)*t + phadd(2) + 2*pi*(h.pitch(1)+df_freq)/h.Settings.fs);end
   % end
   mono=0;
    if isfield(h.Settings,'monostereo')
        if strcmp(h.Settings.monostereo,'mono')
            mono=1;
        end
    end
    
    % if only one channel has so far been defined:
    if df==0 || mono
        h.mwav(chan,:,:) = repmat(h.mwav(chan(1),:,:),length(chan),1); 
    end
    % monaural beats
    if isfield(h.Settings,'monaural')
        if h.Settings.monaural && size(h.mwav,1)==2
            h.mon=h.mwav;
            h.mon(chan(1),:,:) = h.mon(chan(1),:,:) - h.mon(chan(2),:,:);
            h.mon(chan(2),:,:) = h.mon(chan(1),:,:);
        end
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
    if isfield(h.Settings,'attenchan')
        if ismember(chan,h.Settings.attenchan)
            if isfield(h,'vol_atten')
                inten_atten = str2double(get(h.vol_atten,'string'));
            else
                inten_atten = h.Settings.atten; 
            end
            if ~adaptive 
                h.inten_atten = inten_atten;
            elseif adaptive && oddball
                if strcmp(h.Settings.oddballmethod,'intensity')
                    h.inten_atten = [inten_atten, (inten_atten+varlevel)];
                    h.inten_atten = h.inten_atten(h.Seq.signal(tr));
                else
                    h.inten_atten = inten_atten; 
                end
            end
        else
            h.inten_atten = 0; 
        end
    else
        h.inten_atten = 0; 
    end
    % find max rms_sound_dB
    %h.mwav_orig = h.mwav; % non-attenuated version
    temp = reshape(permute(h.mwav,[2,1,3]),length(tsum),[]);
    for i = 1:size(temp,2)
        rms_sound_dB(i) = norm(temp(:,i))/sqrt(length(temp(:,i)));
    end
    rms_sound_dB = max(rms_sound_dB);
    ratio = min(1,(10^(h.inten_atten/20))/rms_sound_dB); % should always be smaller than 1
    h.mwav = ratio*h.mwav;
    if isfield(h,'mon')
        h.mon = ratio*h.mon;
    end
    
    %if size(h.mwav,3)>1
    %    if isfield(h.Settings,'stimdurtype')
    %        if strcmp(h.Settings.stimdurtype,'rand')
    %            h.dur = h.dur(randperm(length(h.dur)));
    %        end
    %    end
    %    % create mask for temporal patterns
    %    mask = [];
    %    for i = 1:size(h.mwav,3)
    %        mask = [mask; i*ones(round(h.dur(i)*h.Settings.fs),1)]; 
    %    end
    %    % splice them together
    %    tempall = [];
    %    tempall2 = [];
    %    for i = 1:size(h.mwav,3)
    %        temp = h.mwav(:,mask==i,i);
    %        temp2 = h.mwav_orig(:,mask==i,i);
    %        if isfield(h,'mon'); tempm = h.mon(:,mask==i,i);end
    %        % apply tapering
    %        if isfield(h.Settings,'Tukey')
    %            if h.Settings.Tukeytype==1
    %                temp = temp.*repmat(tukeywin(size(temp,2),h.Settings.Tukey)',size(temp,1),1);
    %                if isfield(h,'mon'); tempm = tempm.*repmat(tukeywin(size(tempm,2),h.Settings.Tukey)',size(tempm,1),1);end;
    %            end
    %        end
    %        tempall(:,mask==i) = temp;
    %        tempall2(:,mask==i) = temp2;
    %        if isfield(h,'mon'); tempallm(:,mask==i) = tempm;end
    %    end
    %    h.mwav=tempall;
    %    h.mwav_orig=tempall2;
    %    if isfield(h,'mon');h.mon=tempallm;end
    %    
    %end
    if isfield(h.Settings,'Tukey')
        if h.Settings.Tukeytype==2
            h.mwav = h.mwav.*repmat(tukeywin(size(h.mwav,2),h.Settings.Tukey)',size(h.mwav,1),1);
            if isfield(h,'mon');h.mon = h.mon.*repmat(tukeywin(size(h.mon,2),h.Settings.Tukey)',size(h.mon,1),1);end
        end
    end
        
    
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
    if isfield(h,'mon'); 
        h.Seq.stimseq = [h.Seq.stimseq h.mon];
    else
        h.Seq.stimseq = [h.Seq.stimseq h.mwav];
    end;
    
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
                h.Seq.trialend(tr,1) = h.Seq.trialend(tr-1,1)+size(h.mwav,2);
            catch
                h.Seq.trialend(tr,1) = size(h.mwav,2);
            end
            h.Seq.PatternSamples{tr,1} = cumsum(h.dur*h.Settings.fs);
        end 
    end
    
    %figure;plot(h.Seq.stimseq(2,h.Seq.PatternSamples{tr,1}(1)-500:h.Seq.PatternSamples{tr,1}(1)+500));
end
