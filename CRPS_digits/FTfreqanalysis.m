function EEG = FTfreqanalysis(EEG,freqrange,intimes,basetimewin,ncycles)

cfg = [];
cfg.output     = 'pow';
cfg.foi     = freqrange;
cfg.toi          = intimes;
cfg.pad        = 'nextpow2';
cfg.method     = 'mtmconvol';
cfg.taper       = 'hanning';
cfg.t_ftimwin    = ncycles./cfg.foi;  % cycles per time window
%cfg.method     = 'wavelet';                
%cfg.width      = 1;  
cfg.keeptrials = 'yes';
%cfg.precision = 'single';
EEG = ft_freqanalysis(cfg,EEG); 

%% baseline normalisation
cfg=[];
cfg.baseline     = [basetimewin(1) basetimewin(end)];
cfg.baselinetype = 'db';
EEG = ft_freqbaseline(cfg, EEG);

% remove baseline from average (baseline already
% removed from single trials
%baseidx = dsearchn(EEG.time',cfg.baseline')';
%meanVals = repmat(nanmean(inddata(:,baseidx(1):baseidx(2),:), 2), [1 size(inddata,2) 1]);
%inddata = inddata-meanVals;
%meanVals = repmat(nanmean(evdata(:,baseidx(1):baseidx(2),:), 2), [1 size(evdata,2) 1]);
%evdata = evdata-meanVals;