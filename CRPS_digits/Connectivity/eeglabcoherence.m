function eeglabcoherence(basename)

if isunix
    filepath = '/scratch/cb802/Data/CRPS_resting/EEG';
else
    filepath = 'W:\Data\CRPS_resting\EEG';
end

EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);
load chancross;
%freqlist = [0 4; 4 8; 8 13; 13 30; 30 40];
freqlist = [0 4];

chanlocs = EEG.chanlocs;
matrix_max=zeros(size(freqlist,1),EEG.nbchan,EEG.nbchan); % size(freqlist,1) lines ; EEG.nbchan columns ; EEG.nbchan time this table
coh_max = zeros(EEG.nbchan,EEG.nbchan);
cohboot_max = zeros(EEG.nbchan,EEG.nbchan);
matrix_mean=zeros(size(freqlist,1),EEG.nbchan,EEG.nbchan); % size(freqlist,1) lines ; EEG.nbchan columns ; EEG.nbchan time this table
coh_mean = zeros(EEG.nbchan,EEG.nbchan);
cohboot_mean = zeros(EEG.nbchan,EEG.nbchan);
bootrep = 50;


for cc = [1:size(chancross,1)]
    chan1 = chancross{cc,1};
    chan2 = chancross{cc,2};
    chan1i = find(strcmp({chanlocs.labels},chan1));
    chan2i = find(strcmp({chanlocs.labels},chan2));
    [a,cohcc,b,c,freqsout,d,e,f,g,h,cohbootall] = ...
        evalc(['pop_newcrossf( EEG, 1, chan1i, chan2i, [EEG.xmin  EEG.xmax] * 1000, 0 ,' ...
        '''type'', ''wpli'',''freqs'',[0.5 4],''winsize'',EEG.pnts/2,''boottype'',''rand'',''alpha'',0.0002,''padratio'', 1,''plotamp'',''off'',''plotphase'',''off'')']);
    
    wpli_mean(cc,:) = mean(cohcc,2); % Freq x datapoints for each electrode pair
    wpli_boot_mean(cc,:,:) = mean(cohbootall,3); 
    wpli_max(cc,:) = max(cohcc,[],2); % Freq x datapoints for each electrode pair
    wpli_boot_max(cc,:,:) = max(cohbootall,[],3); 
    cc*100/size(chancross,1)
end


for f = 1:size(freqlist,1)
    [M, bstart] = min(abs(freqsout-freqlist(f,1)));
    [M, bend] = min(abs(freqsout-freqlist(f,2)));
    
    coh_mean(:) = 0;
    coh_mean(logical(tril(ones(size(coh_mean)),-1))) = max(wpli_mean(:,bstart:bend),[],2);
    coh_mean = tril(coh_mean,1)+tril(coh_mean,1)';
    matrix_mean(f,:,:) = coh_mean;
    
    
    for r = 1:bootrep
        cohboot_mean(:) = 0;
        cohboot_mean(logical(tril(ones(size(cohboot_mean)),-1))) = max(wpli_boot_mean(:,r,bstart:bend),[],3);
        cohboot_mean = tril(cohboot_mean,1)+tril(cohboot_mean,1)';
        bootmat_mean(f,:,:,r) = cohboot_mean;
    end
    
    coh_max(:) = 0;
    coh_max(logical(tril(ones(size(coh_max)),-1))) = max(wpli_max(:,bstart:bend),[],2);
    coh_max = tril(coh_max,1)+tril(coh_max,1)';
    matrix_max(f,:,:) = coh_max;
    
    
    for r = 1:bootrep
        cohboot_max(:) = 0;
        cohboot_max(logical(tril(ones(size(cohboot_max)),-1))) = max(wpli_boot_max(:,r,bstart:bend),[],3);
        cohboot_max = tril(cohboot_max,1)+tril(cohboot_max,1)';
        bootmat_max(f,:,:,r) = cohboot_max;
    end
end

save(fullfile(filepath, 'eeglab_dwpli', [basename '_eeglab_dwpli.mat']),'matrix_mean','bootmat_mean','matrix_max','bootmat_max','chanlocs');
