function EEG = flipchanEGI(EEG)

load GSN-HydroCel-128-Flipmap.mat

chanlabels = {EEG.chanlocs.labels};

for f = 1:size(flipmap,1)
    chan1idx = strcmp(flipmap{f,1},chanlabels);
    chan2idx = strcmp(flipmap{f,2},chanlabels);
    if sum(chan1idx) == 1 && sum(chan2idx) == 1
        tmpdata = EEG.data(chan1idx,:,:);
        EEG.data(chan1idx,:,:) = EEG.data(chan2idx,:,:);
        EEG.data(chan2idx,:,:) = tmpdata;
    end
end
