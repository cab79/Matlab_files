function EEG = FTfreqanalysis(EEG,S)

if ~isfield(S,'analysistype')
    error('update calling function to include S.analysistype')
end

cfg = [];

if ischar(S.pad)
    if strcmp(S.pad,'min_pad')
        timewin = (S.intimes(end)-S.intimes(1))/1000;
        divi = min(divisors(min(min(S.freqsrange))));
        cfg.pad = divi*ceil(timewin/divi);%-timewin;%'nextpow2';
    else
        cfg.pad = S.pad;
    end
else
    cfg.pad = S.pad;
end

if size(S.freqsrange,1)>1
    cfg.foilim        = [min(S.freqsrange) max(S.freqsrange)];
else
    cfg.foi     = S.freqsrange;
end

if strcmp(S.analysistype,'Freq')
    cfg.output     = 'pow';
    cfg.method     = S.freq_method;
    
elseif strcmp(S.analysistype,'TF')
    cfg.output     = 'pow';
    cfg.method     = S.freq_method;
    cfg.width      = S.freq_param; % frequency-dependent time window
    if cfg.width>0
        cfg.t_ftimwin    = cfg.width./cfg.foi;  % cycles per time window (wavelets only)     
    else
        cfg.t_ftimwin    = min([cfg.pad, 1/min(cfg.foi), 1/min(diff(cfg.foi))])*ones(1,length(cfg.foi)); % works out minimum freq resolution possible
    end
    %cfg.t_ftimwin    = cfg.width./min(cfg.foilim);
    cfg.tapsmofrq    = S.freq_param; % only used for dpss           
     
elseif strcmp(S.analysistype,'Coh')
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
end


cfg.taper       = S.freq_taper;
cfg.toi          = S.intimes/1000;

cfg.keeptrials = 'yes';
%cfg.precision = 'single';
EEG = ft_freqanalysis(cfg,EEG); 

close all
cfg              = [];
cfg.baseline     = 'no';%S.basewin; 
%cfg.baselinetype = 'absolute'; 
cfg.maskstyle    = 'saturation';	
cfg.channel      = 18;
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, EEG);
drawnow