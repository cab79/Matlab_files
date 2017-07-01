
%load chancross; %generated from ft_frequencyanalysis resulting structure EEG.labelcmb
% OR generate chancross from EEG.chanlocs

chans = {EEG.chanlocs.labels};
chancross = repmat(chans,1,length(chans));
chancross(2,:) = reshape(repmat(chans',1,length(chans))',1,length(chans)*length(chans));
chancross=chancross';
allcc=zeros(size(chancross,1),1);
for c = 1:size(chans,2)
    allc=strcmp(chancross(:,2),chans(c));
    allc_ind=find(allc==1);
    allc(allc_ind(1:c)) = zeros(length(allc_ind(1:c)),1);
    allcc=allcc+allc;
end
chancross=chancross(find(allcc==1),:);

%%

no_chan = length(unique(chancross));
chandist = zeros(size(EEG.chanlocs,1),size(EEG.chanlocs,1));

for chan = 1:size(chancross,1)
    chan1 = chancross{chan,1};
    chan2 = chancross{chan,2};
    chan1i = find(strcmp({EEG.chanlocs.labels},chan1));
    chan2i = find(strcmp({EEG.chanlocs.labels},chan2));
    x1 = EEG.chanlocs(chan1i).X;
    y1 = EEG.chanlocs(chan1i).Y;
    z1 = EEG.chanlocs(chan1i).Z;
    x2 = EEG.chanlocs(chan2i).X;
    y2 = EEG.chanlocs(chan2i).Y;
    z2 = EEG.chanlocs(chan2i).Z;
    
    chan1i = find(strcmp({EEG.chanlocs.labels},chan1));
    chan2i = find(strcmp({EEG.chanlocs.labels},chan2));
    
    chandist(chan1i,chan2i) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
    chandist(chan2i,chan1i) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
end

save chandist chandist
save chancross chancross