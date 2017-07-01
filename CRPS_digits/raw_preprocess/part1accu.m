%calculates accuracy for each trial
accudata = [];
for e = 1:length(EEG.event)
    if strcmp(EEG.event(e).type,'STIM')
        thisfnum = EEG.event(e).codes{strcmp('FNUM',EEG.event(e).codes(:,1)),2};
        if e==length(EEG.event)
            thisrnum = thisfnum;
        else
            thisrnum = EEG.event(e+1).codes{strcmp('RNUM',EEG.event(e+1).codes(:,1)),2};
        end
        accu = double(thisfnum == thisrnum);
        if sum(strcmp('ACCU',EEG.event(e).codes(:,1))) == 0
            EEG.event(e).codes = cat(1,EEG.event(e).codes,{'ACCU' accu});
        else
            EEG.event(e).codes{strcmp('ACCU',EEG.event(e).codes(:,1)),2} = accu;
        end
        accudata = cat(1,accudata,[thisfnum accu]);
    end
end
%xlswrite(sprintf('%s_accu.xls',EEG.setname),accudata);
%ALLEEG(CURRENTSET) = EEG;