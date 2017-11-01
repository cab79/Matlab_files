function EEG = FTfreqanalysis(EEG,select,freqrange,intimes,basetimewin,ncycles)

if ~exist('select','var')
    error('update calling function to include select')
end

cfg = [];
if strcmp(select,'Freq')
    cfg.output     = 'pow';
    cfg.method     = 'mtmfft';
elseif strcmp(select,'TF')
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
elseif strcmp(select,'Coh')
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
end

cfg.pad        = 'nextpow2';
if size(freqrange,1)>1
    cfg.foilim        = [min(freqrange) max(freqrange)];
else
    cfg.foi     = freqrange;
end

cfg.toi          = intimes;
cfg.taper       = 'hanning';
cfg.t_ftimwin    = ncycles./cfg.foi;  % cycles per time window
%cfg.method     = 'wavelet';                
%cfg.width      = 1;  
cfg.keeptrials = 'yes';
%cfg.precision = 'single';
EEG = ft_freqanalysis(cfg,EEG); 


if strcmp(select,'TF')
    %% baseline normalisation
    cfg=[];
    cfg.baseline     = [basetimewin(1) basetimewin(end)];
    cfg.baselinetype = 'db';
    EEG = ft_freqbaseline(cfg, EEG);
end

% remove baseline from average (baseline already
% removed from single trials
%baseidx = dsearchn(EEG.time',cfg.baseline')';
%meanVals = repmat(nanmean(inddata(:,baseidx(1):baseidx(2),:), 2), [1 size(inddata,2) 1]);
%inddata = inddata-meanVals;
%meanVals = repmat(nanmean(evdata(:,baseidx(1):baseidx(2),:), 2), [1 size(evdata,2) 1]);
%evdata = evdata-meanVals;