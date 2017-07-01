clear all
loadpaths
loadsubj
grplist = [39 40 41 42];
cond = [1 5];
altcond = [10 6];
sd_correct_cond = [1 2 3 4 5];
altsd_correct_cond = [10 9 8 7 6];
if ~isempty(sd_correct_cond); sdon=1; else sdon=0; end;

subjects = subjlists(grplist);
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        filename = [subj '.set'];
        [pth nme ext] = fileparts(filename);
        TEEG = pop_loadset(filename,filepath);

        dsize = size(TEEG.data);
        %TEEG = pop_eegfilt(TEEG,0,25,0,0);
        TEEG.data = reshape(TEEG.data,dsize(1),dsize(2),dsize(3));

        events = zeros(1,size(TEEG.data,3));
        for i = 1:size(TEEG.data,3)
             if length(TEEG.epoch(i).eventcodes)==1 || length(TEEG.epoch(i).eventcodes)==3
                events(1,i) = TEEG.epoch(i).eventcodes{1,1}{2,2};
             elseif length(TEEG.epoch(i).eventcodes)==2 || length(TEEG.epoch(i).eventcodes)==4
                events(1,i) = TEEG.epoch(i).eventcodes{1,2}{2,2};
             end
        end
        
        if sdon==1
            lscc = length(sd_correct_cond);
            
            for c = 1:lscc
                ceventsL = find(events==sd_correct_cond(c));
                ceventsR = find(events==altsd_correct_cond(c));
                if length(ceventsL) > length(ceventsR)
                    cevents = ceventsL;
                else cevents = ceventsR;
                end
                EEG = pop_select(TEEG, 'trial', cevents);

                for i = 1:size(EEG.data,3)
                     if length(EEG.epoch(i).eventcodes)==1 || length(EEG.epoch(i).eventcodes)==3
                        EEG.epoch(i).eventcodes{1,1}{2,2} = sd_correct_cond(c);
                     elseif length(EEG.epoch(i).eventcodes)==2 || length(EEG.epoch(i).eventcodes)==4
                        EEG.epoch(i).eventcodes{1,2}{2,2} = sd_correct_cond(c);
                     end
                end
                avg = squeeze(mean(EEG.data,3));
                if ~exist('sddat','var'); sddat = zeros(size(avg)); end;
                sddat = cat(3,sddat,avg);
            end
            sd = std(sddat,0,3);
        end

        for c = 1:length(cond)
            ceventsL = find(events==cond(c));
            ceventsR = find(events==altcond(c));
            if length(ceventsL) > length(ceventsR)
                cevents = ceventsL;
            else cevents = ceventsR;
            end
            EEG = pop_select(TEEG, 'trial', cevents);

            for i = 1:size(EEG.data,3)
                 if length(EEG.epoch(i).eventcodes)==1 || length(EEG.epoch(i).eventcodes)==3
                    EEG.epoch(i).eventcodes{1,1}{2,2} = cond(c);
                 elseif length(EEG.epoch(i).eventcodes)==2 || length(EEG.epoch(i).eventcodes)==4
                    EEG.epoch(i).eventcodes{1,2}{2,2} = cond(c);
                 end
            end
            sdmat = reshape(repmat(sd,1,size(EEG.data,3)),size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
            EEG.data = EEG.data./sdmat;
            efilename = [nme '.c' num2str(cond(c)) '.sdnorm' num2str(sdon) ext];
            pop_saveset(EEG, efilename, filepath);
            clear EEG 
        end

        clear TEEG ALLEEG
    end
end

