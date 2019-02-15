function EEGdata = flipchan(EEGdata,chanlocs,varargin)

if nargin>2
    S=varargin{1};
    try
        load(S.path.GSNlocs);
    catch
        load('C:\Data\CORE\eeg\GSN-HydroCel-128-Flipmap.mat');
    end
else
    load('C:\Data\CORE\eeg\GSN-HydroCel-128-Flipmap.mat');
end

chanlabels = {chanlocs.labels};

chanflipped = [];
for f = 1:size(flipmap,1)
    chan1idx = strcmp(flipmap{f,1},chanlabels);
    chan2idx = strcmp(flipmap{f,2},chanlabels);
    if sum(chan1idx) == 1 && sum(chan2idx) == 1
        chanflipped = [chanflipped find(chan1idx) find(chan2idx)];
        tmpdata1 = EEGdata(chan1idx,:,:);
        tmpdata2 = EEGdata(chan2idx,:,:);
        EEGdata(chan1idx,:,:) = tmpdata2;
        EEGdata(chan2idx,:,:) = tmpdata1;
    end
end
disp(['flipped ' num2str(length(chanflipped)) ' channels: ' num2str(sort(chanflipped))]);
