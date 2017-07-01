clear all


%grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\';%Exp1 left v right
grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';%Exp1 left v righ

run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
digit_ranges = {sort(1:5,'descend');6:10};

cd(filepath)
sublist = {};
subjects = subjlists(grplist);
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        accudata = [];
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        [pth nme ext] = fileparts(EEG.filename);
        [pth orig_datfile ext] = fileparts(EEG.filename)

        for e = 1:length(EEG.urevent)
            if strcmp(EEG.urevent(e).type,'STIM')
                thisfnum = EEG.urevent(e).codes{strcmp('FNUM',EEG.urevent(e).codes(:,1)),2};
                r = find(cellfun(@(x)ismember(thisfnum,x),digit_ranges, 'UniformOutput', 1));
                ri = digit_ranges{r};
                nri = 1:2;
                nri(nri==r)=[];
                nr = digit_ranges{nri};
                
                if e==length(EEG.urevent)
                    thisrnum = thisfnum;
                else
                    try
                        thisrnum = EEG.urevent(e+1).codes{strcmp('RNUM',EEG.urevent(e+1).codes(:,1)),2};
                    catch
                        thisrnum = NaN;
                    end
                end
                
                if ~ismember(thisrnum,ri) && thisrnum~=0 && ~isnan(thisrnum)
                    thisrnum=ri(nr==thisrnum);
                end

                if ~isnan(thisrnum)
                    accu = double(thisfnum == thisrnum);
                else
                    accu = NaN;
                end
                
                if sum(strcmp('ACCU',EEG.urevent(e).codes(:,1))) == 0
                    EEG.urevent(e).codes = cat(1,EEG.urevent(e).codes,{'ACCU' accu});
                else
                    EEG.urevent(e).codes{strcmp('ACCU',EEG.urevent(e).codes(:,1)),2} = accu;
                end
                accudata = cat(1,accudata,[thisfnum thisrnum accu]);
            end
        end
        
        for e = 1:length(EEG.event)
            if strcmp(EEG.event(e).type,'STIM')
                if sum(strcmp('ACCU',EEG.event(e).codes(:,1))) == 0
                    EEG.event(e).codes = cat(1,EEG.event(e).codes,{'ACCU' EEG.urevent(EEG.event(e).urevent).codes{strcmp('ACCU',EEG.urevent(EEG.event(e).urevent).codes(:,1)),2}});
                else
                    EEG.event(e).codes{strcmp('ACCU',EEG.event(e).codes(:,1)),2} = EEG.urevent(EEG.event(e).urevent).codes{strcmp('ACCU',EEG.urevent(EEG.event(e).urevent).codes(:,1)),2};
                end
            end
        end
        for e = 1:length(EEG.epoch)
            if sum(strcmp(EEG.epoch(e).eventtype(1,:),'STIM'))>0
                ind = find(strcmp(EEG.epoch(e).eventtype(1,:),'STIM'));
                for inds = 1:length(ind)
                    findaccu = strcmp('ACCU',EEG.epoch(e).eventcodes{1,ind(inds)}(:,1));
                    if sum(findaccu) > 0
                        EEG.epoch(e).eventcodes{1,ind(inds)}(find(findaccu),:) = [];
                    end
                    EEG.epoch(e).eventcodes{1,ind(inds)} = cat(1,EEG.epoch(e).eventcodes{1,ind(inds)},{'ACCU' EEG.urevent(EEG.epoch(e).eventinit_index{1,ind(inds)}).codes{strcmp('ACCU',EEG.urevent(EEG.epoch(e).eventinit_index{1,ind(inds)}).codes(:,1)),2}});
                end
            end
        end
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',EEG.filename,'filepath',filepath); 
        %xlswrite(sprintf('%s_accu.xls',EEG.setname),accudata);
    end
end