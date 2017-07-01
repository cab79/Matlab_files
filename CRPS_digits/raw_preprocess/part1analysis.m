
selectaccu = [0 1];

fnumepochs = false(1,length(EEG.epoch));
accuepochs = false(1,length(EEG.epoch));
trials=cell(1,length(EEG.epoch));
for ep = 1:length(EEG.epoch)
    stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
    if ~isempty(stimevidx)
        stimcodes = EEG.epoch(ep).eventcodes{stimevidx};
        
        thisfnum = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
 %       thisaccu = stimcodes{strcmp('ACCU',stimcodes(:,1)),2};
        
        if sum(thisfnum == selectfnum) == 1
            fnumepochs(ep) = 1;
        end
        
 %       if sum(thisaccu == selectaccu) == 1
 %           accuepochs(ep) = 1;
 %       end
        
        selectepochs = fnumepochs;%(fnumepochs & accuepochs);
    end
    trials{1,ep} = num2str(thisfnum); 
end

EEG = pop_select(EEG,'trial',find(selectepochs));
EEG.conditionlabel = trials(fnumepochs);