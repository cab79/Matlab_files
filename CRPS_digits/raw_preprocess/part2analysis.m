%selectfnum =1:5;

selectedepochs = false(1,length(EEG.epoch));
trials=cell(1,length(EEG.epoch));
for ep = 1:length(EEG.epoch)
    %if size(EEG.epoch(ep).eventcodes,2)==2; continue; end
    thisfnum = EEG.epoch(ep).eventcodes{find(strcmp('FNUM',EEG.epoch(ep).eventcodes(:,1))),2};
    if sum(thisfnum == selectfnum) == 1
        selectedepochs(ep) = 1;
    end 
    
    trials{1,ep} = num2str(thisfnum); 
end

EEG = pop_select(EEG,'trial',find(selectedepochs));

EEG.conditionlabel = trials(selectedepochs);