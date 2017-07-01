function ftcoherence(basename)
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
cohboot = zeros(EEG.nbchan,EEG.nbchan);
bootrep = 50;


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

EEGf = ft_freqanalysis(cfg,EEG); % returns power (auto-spectra) and shared activity (cross-spectra) between all channels
wpli = ft_connectivity_wpli(EEGf.crsspctrm,'debias',true,'dojack',false);

for f = 1:size(freqlist,1)
    fprintf('f %d',f);
    [M, bstart] = min(abs(EEGf.freq-freqlist(f,1)));
    [M, bend] = min(abs(EEGf.freq-freqlist(f,2)));
    
    coh(:) = 0;
    coh(logical(tril(ones(size(coh)),-1))) = max(wpli(:,bstart:bend),[],2);
    coh = tril(coh,1)+tril(coh,1)';
    matrix(f,:,:) = coh;
end
clear EEGf coh wpli

for nboot = 1:bootrep
    fprintf('nboot %d',nboot);
    for tr = 1:length(EEG.trial)
        for ele = 1:size(EEG.trial{tr},1);
            trial = EEG.trial{tr};
            surrEEG=[];
            surrEEG = [phaseran(squeeze(trial(tr,:)'),1); 0]';
            %surrEEG = [phaseran(EEG.trial{1}',1); zeros(1, length(EEG.label))]';
            trial(tr,:) = surrEEG;
        end
        EEG.trial{tr} = trial;
    end
    
    
    EEGf = ft_freqanalysis(cfg,EEG); % returns power (auto-spectra) and shared activity (cross-spectra) between all channels
    wpli_boot = ft_connectivity_wpli(EEGf.crsspctrm,'debias',true,'dojack',false);
    
    for f = 1:size(freqlist,1)
        [M, bstart] = min(abs(EEGf.freq-freqlist(f,1)));
        [M, bend] = min(abs(EEGf.freq-freqlist(f,2)));

        cohboot(:) = 0;
        cohboot(logical(tril(ones(size(cohboot)),-1))) = max(wpli_boot(:,bstart:bend),[],2);
        cohboot = tril(cohboot,1)+tril(cohboot,1)';
        bootmat(f,:,:,nboot) = cohboot;
    end
    clear EEGf wpli_boot
end


save(fullfile(filepath, 'ft_dwpli', [basename '_ftdwpli.mat']),'matrix','bootmat','chanlocs');
