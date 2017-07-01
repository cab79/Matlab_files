function EEG = sortchan(EEG)

chanlabels = {EEG.chanlocs.labels};
channum = zeros(size(chanlabels));

for c = 1:length(chanlabels)
    cnum = str2double(chanlabels{c}(2:end));
    if ~isempty(cnum)
        channum(c) = cnum;
    end
end

[channum,sortidx] = sort(channum);
sortidx = [sortidx(channum ~= 0) sortidx(channum == 0)];
EEG.chanlocs = EEG.chanlocs(sortidx);
EEG.data = EEG.data(sortidx,:,:);