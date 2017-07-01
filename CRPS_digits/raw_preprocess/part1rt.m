edgethresh = 350;

emgidx = find(strcmp('EMG',{EEG.chanlocs.labels}));
for e = 1:length(EEG.event)
    switch(EEG.event(e).type)
        case 'STIM'
            EEG.data(emgidx,EEG.event(e).latency:EEG.event(e).latency+round(EEG.srate/2)) = 0;
        case 'RESP'
            EEG.data(emgidx,EEG.event(e).latency-round(EEG.srate/2):EEG.event(e).latency) = 0;
    end
end

EEG = pop_chanevent(EEG, emgidx,'oper',sprintf('X > %d | X < %d', edgethresh, -edgethresh),...
    'edge','leading','edgelen',EEG.srate,'delchan','off','delevent','off','nbtype',1,'typename','RT');

rtdata = [];
for e = 2:length(EEG.event)
    if strcmp(EEG.event(e).type,'RT') && strcmp(EEG.event(e-1).type,'STIM')
        trialrt = (( EEG.event(e).latency - EEG.event(e-1).latency ) * 1000)/ EEG.srate;
        if sum(strcmp('RT',EEG.event(e-1).codes(:,1))) == 0
            EEG.event(e-1).codes = cat(1,EEG.event(e-1).codes,{'RT',trialrt});
        else
            EEG.event(e-1).codes{strcmp('RT',EEG.event(e-1).codes(:,1)),2} = trialrt;
        end
        
        thisfnum = EEG.event(e-1).codes{strcmp('FNUM',EEG.event(e-1).codes(:,1)),2};
        rtdata = cat(1,rtdata,[thisfnum trialrt]);
    end
end

xlswrite(sprintf('%s_rt.xls',EEG.setname),rtdata);

ALLEEG(CURRENTSET) = EEG;
eeglab redraw