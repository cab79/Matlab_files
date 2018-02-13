function IAF=autoIAF

% INSTRUCTIONS:
% put data to be analysed in a folder and put the folder name below.
% The script will analyse EVERY file in the folder, so only have files in
% there that require analysis (otherwise move them)

% SET SOME OPTIONS
S.filepath = 'C:\Data\NTIP\Preprocessed';
S.file_ext = '.set'; % set of vhdr
S.preprocess = 0; 
S.filterset = [6 14]; % FILTER AND PLOT SETTINGS
S.freqwin = [8 12]; % window for peak detection
S.removechan = [32]; % ECG
S.timebin = [0 0.5]; % for epoching
S.epochrejthresh = 1000; % data with very large deviations will be removed

% RUN
dbstop if error
close all
IAF = runautoIAF(S);
end

function IAF = runautoIAF(S)

% CONVERT
filepath = S.filepath;
fnames = dir(fullfile(filepath,['*' S.file_ext]));
IAF = nan(length(fnames),1);
for fn = 1:length(fnames)
    [~,filename,~] = fileparts(fnames(fn).name);
    if strcmp(S.file_ext,'.vhdr')
        dataimport_bv(filepath,[filename '.vhdr'],1,0,filepath)
        setfile = [filename '_orig.set'];
    else
        setfile = [filename '.set'];
    end

    % LOAD EEGLAB DATA
    EEG = pop_loadset('filename',setfile,'filepath',filepath);

    if S.preprocess
    
        % remove chans
        EEG = pop_select(EEG,'nochannel',S.removechan);

        % FILTER
        if S.filterset(1)>0; EEG = pop_eegfiltnew( EEG, S.filterset(1), 0, [], 0);end
        if S.filterset(2)>0; EEG = pop_eegfiltnew( EEG, 0, S.filterset(2), [], 0);end

        % EPOCH
        %add 2sec epochs by adding markers every second and epochs after that
        %add the marker ('M') 
        Sr = EEG.srate; % sampling rate of data
        Ndp = Sr*(S.timebin(2)-S.timebin(1));% number of data points per epoch
        Tdp = size(EEG.data,2);% total number of data points in file
        Mep = floor(Tdp/Ndp);% max possible number of epochs
        for i = 1:Mep;
            EEG.event(1,i).type = 'M';
            EEG.event(1,i).latency = (i-1)*Ndp+1;
            EEG.event(1,i).urevent = 'M';
        end

        %create 2sec epochs
        EEG = pop_epoch( EEG, {  'M'  }, S.timebin);

        % LINEAR DETREND
        %for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;

        % REMOVE BASELINE
        %EEG = pop_rmbase( EEG, [timebin(1)*100    0]);

        %REJECT CHANNELS
        [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'measure','kurt','freqrange',S.filterset);

        %REJECT EPOCHS
        EEG = pop_autorej(EEG, 'nogui','on','threshold',S.epochrejthresh,'startprob',12);
    end

    %FREQ SPECTRUM
    y = reshape(EEG.data,EEG.nbchan,{})';
    Ny = length(y);
    w = hanning(Ny); % Window data
    nfft = Ny;
    psd=[];
    for i= 1:size(y,2)
        [psd(:,i),freq] = pwelch(y(:,i),w,0,nfft,EEG.srate);
    end
    %psd = fft(y,nfft)/Ny;
    psd = mean(psd,2);
    freq = EEG.srate/2*linspace(0,1,nfft/2+1);
    selectfreq = dsearchn(freq',S.filterset');
    freq=freq(selectfreq(1):selectfreq(2));
    psd =psd(selectfreq(1):selectfreq(2));

    %yy1 = smooth(freq,psd,0.15,'loess');
    yy2 = smooth(freq,psd,0.15,'rloess');
    figure; plot(freq,10*log10(psd),'Color', [0.8 0.8 0.8]);
    hold on
    %plot(freq,10*log10(yy1),'r');
    %hold on
    plot(freq,10*log10(yy2),'b');
    hold on

    % find IAF
    freqwindow = dsearchn(freq',S.freqwin');
    [PKS,LOCS]=max(yy2(freqwindow(1):freqwindow(2)));
    IAF(fn,1)=freq(freqwindow(1)+LOCS-1);
    scatter(IAF(fn,1),10*log10(PKS),'m');
    title(filename)
    xlabel('Frequency (Hz)')
    ylabel('Power Spectral Density')
end

end
