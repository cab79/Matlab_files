function EEG = FTfreqanalysis(EEG,S)

if ~isfield(S.(S.func).select,'datatype')
    error('update calling function to include S.select.datatype')
end

cfg = [];

if ischar(S.(S.func).freq.pad)
    if strcmp(S.(S.func).freq.pad,'min_pad')
        timewin = (S.(S.func).intimes(end)-S.(S.func).intimes(1))/1000;
        divi = min(divisors(min(min(S.(S.func).select.freq))));
        cfg.pad = divi*ceil(timewin/divi);%-timewin;%'nextpow2';
    else
        cfg.pad = S.(S.func).freq.pad;
    end
else
    cfg.pad = S.(S.func).freq.pad;
end

if size(S.(S.func).select.freq,1)>1
    cfg.foilim        = [min(S.(S.func).select.freq) max(S.(S.func).select.freq)];
else
    cfg.foi     = S.(S.func).select.freq;
end

if strcmp(S.(S.func).select.datatype,'Freq')
    cfg.output     = 'pow';
    cfg.method     = S.(S.func).freq.method;
    
elseif strcmp(S.(S.func).select.datatype,'TF')
    cfg.output     = 'pow';
    cfg.method     = S.(S.func).freq.method;
    cfg.width      = S.(S.func).freq.param; % frequency-dependent time window
    if cfg.width>0
        cfg.t_ftimwin    = cfg.width./cfg.foi;  % cycles per time window (wavelets only)     
    else
        cfg.t_ftimwin    = min([cfg.pad, 1/min(cfg.foi), 1/min(diff(cfg.foi))])*ones(1,length(cfg.foi)); % works out minimum freq resolution possible
    end
    %cfg.t_ftimwin    = cfg.width./min(cfg.foilim);
    cfg.tapsmofrq    = S.(S.func).freq.param; % only used for dpss           
     
elseif strcmp(S.(S.func).select.datatype,'Coh')
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
end


cfg.taper       = S.(S.func).freq.taper;
cfg.toi          = S.(S.func).intimes/1000;

cfg.keeptrials = 'yes';
%cfg.precision = 'single';
EEG = ft_freqanalysis(cfg,EEG); 

if strcmp(S.(S.func).select.datatype,'TF')
    close all
    cfg              = [];
    cfg.baseline     = 'no';%S.epoch.basewin; 
    %cfg.baselinetype = 'absolute'; 
    cfg.maskstyle    = 'saturation';	
    cfg.channel      = 18;
    cfg.interactive  = 'no';
    figure
    ft_singleplotTFR(cfg, EEG);
end
drawnow