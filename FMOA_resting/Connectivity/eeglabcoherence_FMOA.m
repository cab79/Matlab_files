function eeglabcoherence_FMOA(filepath, basename,basenameb,EEG,conntype)

load(fullfile(filepath,'EEG',[basename '.mat']));
load(fullfile(filepath,'EEG',[basenameb '.mat']));

% exclude bad trials and electrodes
data = total_data_ICA([1:2 4:30 33:64],:,find(events_mat(:,3)==1));
events_mat = events_mat(find(events_mat(:,3)==1),:);
clear total_data_ICA;

% select x number of trials
ntrials=15;
eopen = find(events_mat(:,2)==1);
eopen = eopen(1:min(ntrials,length(eopen)));
eclosed = find(events_mat(:,2)==2);
eclosed = eclosed(1:min(ntrials,length(eclosed)));
data =cat(3, data(:,:,eopen),data(:,:,eclosed));

EEG.chanlocs = EEG.chanlocs([1:2 4:30 33:64]);
[nele npnts nep] = size(data);
data = reshape(data,nele,npnts*nep);
load(fullfile(filepath,'EEG','chancross.mat'));
EEG.srate = 250;
EEG.nbchan = nele;
freqlist = [0 4; 4 8; 8 13; 13 30; 30 40; 60 100];
%freqlist = [0 4];


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
    chan1i = find(strcmp({EEG.chanlocs.labels},chan1));
    chan2i = find(strcmp({EEG.chanlocs.labels},chan2));
    %[a,cohcc,b,c,freqsout,d,e,f,g,h,cohbootall] = ...
    %    evalc(['pop_newcrossf( EEG, 1, chan1i, chan2i, [EEG.xmin  EEG.xmax] * 1000, 0 ,' ...
    %    '''type'', ''wpli'',''winsize'',EEG.pnts/2,''boottype'',''rand'',''alpha'',0.0002,''padratio'', 1,''plotamp'',''off'',''plotphase'',''off'')']);
    
    chan1ts = data(chan1i,:);
    chan2ts = data(chan2i,:);
    [a, cohcc,b,c,freqsout,d,e,f,g,h,cohbootall] = ...
        evalc(['newcrossf( chan1ts, chan2ts, npnts, [0 npnts]*(1000/EEG.srate), EEG.srate, 0 ,' ...
        '''type'', ''wpli'',''winsize'',npnts/2,''freqs'',[0.5 100],''boottype'',''rand'',''alpha'',0.0002,''padratio'', 1,''plotamp'',''off'',''plotphase'',''off'')']);
    %[R,mbase,timesout,freqs,Rbootout,Rangle, coherresout, alltfX, alltfY,Rbootall] = newcrossf(chan1i, chan2i, npnts, [0 npnts]*(1000/EEG.srate), EEG.srate, 0, varargin)
    
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
chanlocs = EEG.chanlocs;
save(fullfile(filepath, 'EEG', conntype, [basename '_' conntype '.mat']),'matrix_mean','bootmat_mean','matrix_max','bootmat_max','chanlocs');
