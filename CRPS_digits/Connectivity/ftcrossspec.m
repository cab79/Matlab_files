function ftcrossspec(basename)
% matrix with all coherence number for a frequency interval : freq1:freq2

if isunix
    filepath = '/scratch/cb802/Data/CRPS_resting/EEG';
else
    filepath = 'W:\Data\CRPS_resting\EEG';
end

%     if exist([filepath basename 'wplifdr.mat'],'file')
%         fprintf('%s exists. Skipping...\n',[basename 'wplifdr.mat']);
%         continue;
%     end

EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);
%load(fullfile(filepath, 'spectra', [basename '_spectra.mat']),'freqlist');
freqlist = [0 4; 4 8; 8 13; 13 30; 30 40];

chanlocs = EEG.chanlocs;
matrix=zeros(size(freqlist,1),EEG.nbchan,EEG.nbchan); % size(freqlist,1) lines ; EEG.nbchan columns ; EEG.nbchan time this table
coh = zeros(EEG.nbchan,EEG.nbchan);

EEG = convertoft(EEG);
cfg = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim        = [0.1 50];
cfg.taper = 'hanning';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 0.5;
cfg.keeptrials = 'yes';
cfg.precision = 'single';

EEG = ft_freqanalysis(cfg,EEG); % returns power (auto-spectra) and shared activity (cross-spectra) between all channels
cros = abs(EEG.crsspctrm);
cros = squeeze(mean(cros,1));

for f = 1:size(freqlist,1)
    [M, bstart] = min(abs(EEG.freq-freqlist(f,1)));
    [M, bend] = min(abs(EEG.freq-freqlist(f,2)));
    
    coh(:) = 0;
    coh(logical(tril(ones(size(coh)),-1))) = max(cros(:,bstart:bend),[],2);
    coh = tril(coh,1)+tril(coh,1)';
    matrix(f,:,:) = coh;
end

save(fullfile(filepath, 'ft_crossspec', [basename '_ftcrossspec.mat']),'matrix','chanlocs');
