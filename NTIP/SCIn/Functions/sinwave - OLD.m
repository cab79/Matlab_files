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

% set default to sin wave
if ~isfield(h.Settings,'wavetype')
    h.Settings.wavetype = 'sin';
elseif isempty(h.Settings.wavetype)
    h.Settings.wavetype = 'sin';
end
% set default to phase alignment
if ~isfield(h.Settings,'alignphase')
    h.alignphase = 1;
elseif isempty(h.Settings.alignphase)
    h.alignphase = 1;
else
    h.alignphase = 0;
end

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

% if THRESHOLD
threshold=0;
if isfield(h.Settings,'threshold')
    if ~isempty(h.Settings.threshold)
        threshold=1;
    end
end

% if ODDBALL
oddball=0;
if isfield(h.Settings,'oddballmethod')
    if ~isempty(h.Settings.oddballmethod)
        oddball=1;
    end
end

% find trial(s) for which to create wave
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
    h.freq = h.Settings.f0;
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
                    if strcmp(conditionmethod,'pitch') || strcmp(conditionmethod,'freq')
                        if iscell(h.Settings.conditionvalue)
                            h.freq = h.Settings.conditionvalue{h.Seq.signal(tr),i};
                        else
                            h.freq = h.Settings.conditionvalue(i,h.Seq.signal(tr));
                        end
                    end
                    if strcmp(conditionmethod,'intensity')
                        if iscell(h.Settings.conditionvalue)
                            h.inten = h.Settings.conditionvalue{h.Seq.signal(tr),i};
                        else
                            h.inten = h.Settings.conditionvalue(i,h.Seq.signal(tr));
                        end
                    end
                    if strcmp(conditionmethod,'phase')
                        if iscell(h.Settings.conditionvalue)
                            h.alignphase = h.Settings.conditionvalue{h.Seq.signal(tr),i};
                        else
                            h.alignphase = h.Settings.conditionvalue(i,h.Seq.signal(tr));
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
        if ~adaptive && ~threshold
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
            elseif strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                h.freq = oddval;
            elseif strcmp(h.Settings.oddballmethod,'duration')
                h.dur = oddval;
            end
        elseif adaptive || threshold
            if isfield(h,'s')
                varlevel = h.s.a(h.Seq.adapttype(h.i)).StimulusLevel;
            else
                if adaptive
                    varlevel = h.Settings.adaptive.startinglevel;
                else
                    varlevel = h.Settings.threshold.startinglevel;
                end  
            end
            if strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                h.freq = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+varlevel)]; % create new pitch pair
                h.freq = h.freq(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddballmethod,'intensity')
                h.inten = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+varlevel)]; % create new pitch pair
                h.inten = h.inten(h.Seq.signal(tr));
            elseif strcmp(h.Settings.oddballmethod,'duration') && (strcmp(h.Settings.patternmethod,'pitch') || strcmp(h.Settings.patternmethod,'freq'))
                if iscell(h.Settings.oddballvalue)
                    h.dur = h.Settings.oddballvalue{h.Seq.signal(tr),:};
                else
                    h.dur = h.Settings.oddballvalue(h.Seq.signal(tr),:);
                end
                h.freq = [h.Settings.patternvalue(1), (h.Settings.patternvalue(1)+varlevel)]; % create new pitch pair
            end
        end
    end
       
    %apply pitch pattern?
    freqpattern=0;
    if isfield(h.Settings,'patternmethod')
        if strcmp(h.Settings.patternmethod,'pitch') || strcmp(h.Settings.patternmethod,'freq') % pitch changes
            freqpattern=1;
            if ~((adaptive || threshold) && (strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq'))) && ~(~isempty(strcmp(h.Settings.conditionmethod,'pitch')) || ~isempty(strcmp(h.Settings.conditionmethod,'freq'))) % pitch already defined above in this case
                if isnumeric(h.Settings.patternvalue)
                    h.freq = h.Settings.patternvalue;
                elseif iscell(h.Settings.patternvalue)
                    nDur = length(h.dur);
                    nPit = cellfun(@length,h.Settings.patternvalue);
                    h.freq = h.Settings.patternvalue{nPit==nDur};
                end
            end
        end
    end
    % apply response probe?
    if isfield(h.Settings,'RPmethod')
        if strcmp(h.Settings.RPmethod,'pitch') || strcmp(h.Settings.RPmethod,'freq')
            if h.Seq.RP(tr)==1
                freqpattern=1;
                h.freq = h.Settings.RPvalue;
                h.dur = h.Settings.RPdur;
            end
        end
    end
    %apply intensity pattern?
    intenpattern=0;
    if isfield(h.Settings,'patternmethod')
        if strcmp(h.Settings.patternmethod,'intensity') % intensity changes
            intenpattern=1;
            if ~any(strcmp(h.Settings.conditionmethod,'intensity')) % then already defined
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
    
    if isfield(h.Settings,'df') && length(chan)==2 % if freq is different in two channels
        df=1;
        if isfield(h,'entrainfreq')
            try
                df_freq = str2double(get(h.entrainfreq,'string'));
            catch
                df_freq = str2double(h.entrainfreq);
            end
            if df_freq==0 || isnan(df_freq)
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
            elseif strcmp(h.Settings.patternmethod,'pitch') || strcmp(h.Settings.patternmethod,'freq')
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
                phadd = zeros(length(chan),1);
            end
        else
            phadd = zeros(length(chan),1);
        end
        
        % initialise
        mwav{i}=zeros(h.Settings.nrchannels,length(t{i}));
        temp_sin{i}=zeros(h.Settings.nrchannels,length(t{i}));
        
        % pitch/inten-specific 
        if freqpattern
            %h.freq = h.freq(rs); % randomise?
            if length(h.dur)~=length(h.freq)
                error('num column of stimdur must equate to number of frequencies');
            end
            usefreq = h.freq(i);
            useinten = h.inten(1);
        elseif intenpattern
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
            mwav{i}(chan(1),:) = useinten *sin(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs);
            if df; mwav{i}(chan(2),:) = useinten *sin(2*pi*(usefreq+df_freq)*t{i} + phadd(2) + 2*pi*(usefreq+df_freq)/h.Settings.fs);end
            temp_sin{i} = mwav{i};
        elseif strcmp(h.Settings.wavetype,'square')
            if isfield(h,'dutycycle')
                if h.dutycycle>0
                    useinten = str2double(h.dutycycle);
                end
            end
            mwav{i}(chan(1),:) = square(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs,useinten);
            temp_sin{i}(chan(1),:) = sin(2*pi*(usefreq)*t{i} + phadd(1) + 2*pi*usefreq/h.Settings.fs); % for phase estimation
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
        h.iphase(:,tr,:,i) = [asin(temp_sin{i}(:,end-1)), temp_sin{i}(:,end)-temp_sin{i}(:,end-1)];% mwav_orig(:,end-1)-mwav_orig(:,end)];
        
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
    if isfield(h.Settings,'attenchan')
        if any(ismember(chan,h.Settings.attenchan))
            if isfield(h,'vol_atten')
                try
                    inten_atten = str2double(get(h.vol_atten,'string'));
                catch
                    inten_atten = str2double(h.vol_atten);
                end
            else
                inten_atten = h.Settings.atten; 
            end
            if ~adaptive && ~threshold
                h.inten_atten = inten_atten;
            elseif (adaptive || threshold) && oddball
                if strcmp(h.Settings.oddballmethod,'intensity')
                    h.inten_atten = [inten_atten, (inten_atten+varlevel)];
                    h.inten_atten = h.inten_atten(h.Seq.signal(tr));
                else
                    h.inten_atten = inten_atten+varlevel;
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
    
    %if isfield(h,'mon')
    %    h.mon = ratio*h.mon;
    %end
    
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
                h.Seq.trialend(tr,1) = h.Seq.trialend(tr-1,1)+size(h.mwav,2);
            catch
                h.Seq.trialend(tr,1) = size(h.mwav,2);
            end
            h.Seq.PatternSamples{tr,1} = cumsum(h.dur*h.Settings.fs);
        end 
    end
    
    %figure;plot(h.Seq.stimseq(2,h.Seq.PatternSamples{tr,1}(1)-500:h.Seq.PatternSamples{tr,1}(1)+500));
end
