clear all
close all

%---settings---%
grplist = [33 34 31 32]; sublist_side = {'A','U','A','U'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';typehand = '_AU';%Exp2
%grplist = [39:42]; sublist_side = {'A','U','A','U'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\'; typehand = '_AU';%Exp1 left v right
%grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\'; typehand = '_LR'; %Exp1 left v righ
%grplist = [43:46]; sublist_side = {'A','U','A','U'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\'; typehand = '_AU'; %Exp1 left v righ

%select = 'TF'; TFmethod = 'FT'; % or, 'EL'
select = 'ERP'; TFmethod = '';
flip=1;
singtrial=0;
eventtypes = {'FNUM','CNUM','ACCU'; [1:5], [0 1], [0 1]}; use_etype = [1 2];
%posnegpeak = [2 1 1 2 2 2 1 1 1 1 1 1 1 1];
numERPcomp=4;
rmCommCond=0;
applyCSD=0; 
if strcmp(select,'ERP') && applyCSD; load C:\Data\Matlab\Matlab_files\CRPS_digits\CSDmontage_92.mat;end
%--------------%

savenametype='';
for i = use_etype
    savenametype = [savenametype '_' eventtypes{1,i}];
end

run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath)
ele_left = [];
ele_right = [];
%ele_left = [29 23 24 28 30 33 34];
%ele_right = [78 66 70 77 79 83 84];

aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1;
    1 2 2 1 2 2 2 1 2 1 2 2 2];
%aff_side = [1 1 1 1 1 1 1 1 1 1 1 1 1];

%cond_nme = {'little','ring','middle','index','thumb'};
hand_nme = {'L','R'};
handaff_nme = {'A','U'};
el = 92;

subjects = subjlists(grplist);
DATs = subjects;

if strcmp(select,'ERP')
    dp = 250;
    nf = 1;
    times = -0.2:0.004:0.796;
elseif strcmp(select,'TF')
    dp = 181;
    freqrange = [4:2:40];
    nf = length(freqrange);
    intimes = -0.2:0.004:0.796;
    basebin = [-0.2 0];
end
    
% create DATs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        chanlocs=EEG.chanlocs;
        dsize = size(EEG.data);
        
        
        %for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
        
        if flip==1 && (~isempty(strfind(subj, 'right')) || ~isempty(strfind(subj, 'Right'))); 
            EEG = flipchan(EEG); 
        end
        
        %EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        condevents = EEG2condevents(EEG,'STIM', eventtypes,use_etype);
        no_cond = length(unique(condevents));
        
        % data cleaning and ERP
        %EEGrep=EEG; % most repeatable over trials
        %EEGrep.data = DSSonERP(EEGrep.data,[2 1 3],condevents,[4],[],'off'); % maximise repeatability
        %pause
        %close all
        DATs{s,1}{s2,1} = cell(no_cond,1);
        SINGDAT = cell(1);
        for i = unique(condevents)
            eegdata = [];
            eegdatasing = [];
            inddata = [];
            inddatasing = [];
            evdata = [];
            evdatasing = [];
            itcdata = [];
            condidx = find(condevents==i);
            if length(condidx)>1
               evdatasing = EEG.data(:,:,condidx);
                if strcmp(select,'ERP') && applyCSD; 
                    [nchan,nsamp,ntrial] = size(evdatasing);
                    evdatasing = reshape(evdatasing,nchan,nsamp*ntrial);
                    evdatasing=CSD(evdatasing,G,H); 
                    evdatasing = reshape(evdatasing,nchan,nsamp,ntrial);
                end
                %evdatasing = DSSonERP(evdatasing,[2 1 3],[],[numERPcomp],[],'off'); % maximise repeatability 
                %trialdata = EEG.data(:,:,condidx);
                evdata = double(mean(evdatasing,3));
                %EEG.data(:,:,condidx) = trialdata;
                if strcmp(select,'TF')
                    trialdata= EEG.data(:,:,condidx);
                    trialdata2= trialdata;
                    for t = 1:size(trialdata,3)
                        trialdata2(:,:,t) = trialdata(:,:,t)-evdata; % subtract out ERP for induced activity
                    end
                    eegdata2 = [];
                    eegdatasing2 = [];
                    if strcmp(TFmethod,'EL')
                        for e = 1:size(EEG.data,1)
                            datvec = reshape(squeeze(trialdata(e,:,:)),1,size(trialdata,2)*size(trialdata,3));
                            datvec2 = reshape(squeeze(trialdata2(e,:,:)),1,size(trialdata2,2)*size(trialdata2,3));
                            %[eegdata(:,:,e),x,y,times,freqs] = timef(datvec,size(EEG.data,2),[-200 800],EEG.srate,1,'detret','on','winsize',60,'plotersp','off','plotitc','off','plotphase','off');  
                            [eegdata(:,:,e),eegdatasing(:,:,:,e),itcdata(:,:,e),~,times,freqs,~,~] = newtimef_singletrial(datvec,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                            [eegdata2(:,:,e),eegdatasing2(:,:,:,e),~,~,times,freqs,~,~] = newtimef_singletrial(datvec2,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                            %[eegdata(:,:,e),eegdatasing(:,:,:,e),~,~,times,freqs] = timef_singletrial(datvec,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'detrend','off','plotersp','off','plotitc','off','plotphase','off');  
                        end
                        inddata = permute(eegdata2,[3 2 1]);
                        evdata = permute(eegdata-eegdata2,[3 2 1]);
                        itcdata = permute(itcdata,[3 2 1]);
                        inddatasing = permute(eegdatasing2,[4 2 3 1]);
                        evdatasing = permute(eegdatasing-eegdatasing2,[4 2 3 1]);
                    elseif strcmp(TFmethod,'FT')
                        FTEEG = pop_select(EEG,'trial',condidx); % previously: selectevent 'event'
                        FTEEG2 = FTEEG;
                        FTEEG2.data = trialdata2;
                        FTEEG2=convertoft(FTEEG2);
                        FTEEG=convertoft(FTEEG);
                        
                        % freq analysis
                        ncycles=1;
                        basetimewin     = [-0.2 0];
                        FTEEG = FTfreqanalysis(FTEEG,freqrange,intimes,basetimewin,ncycles);
                        FTEEG2 = FTfreqanalysis(FTEEG2,freqrange,intimes,basetimewin,ncycles);
                        
                        inddatasing = permute(FTEEG2.powspctrm,[2 4 1 3]);
                        evdatasing = permute(FTEEG.powspctrm,[2 4 1 3])-inddatasing;
                        inddata = squeeze(nanmean(inddatasing,3));
                        evdata = squeeze(nanmean(evdatasing,3));
                        times = FTEEG.time;
                        freqs = freqrange;
                        
                        % remove baseline from average (baseline already
                        % removed from single trials
                        baseidx = dsearchn(times',basebin')';
                        meanVals = repmat(nanmean(inddata(:,baseidx(1):baseidx(2),:), 2), [1 size(inddata,2) 1]);
                        inddata = inddata-meanVals;
                        meanVals = repmat(nanmean(evdata(:,baseidx(1):baseidx(2),:), 2), [1 size(evdata,2) 1]);
                        evdata = evdata-meanVals;
                    end
                end
            end
            
            SINGDAT{i,1} = evdatasing;
            SINGDAT{i,2} = inddatasing;
            SINGDAT{i,3} = length(condidx);
            
            DATs{s,1}{s2,1}{i,1} = evdata;
            DATs{s,1}{s2,1}{i,2} = inddata;
            DATs{s,1}{s2,1}{i,3} = itcdata;
            DATs{s,1}{s2,1}{i,4} = length(condidx);
            %figure
            %plot(EEG.times',evdata);
            %pause
            %close all

        end
        %-- remove maximally similar responses between conditioss --%
        if rmCommCond
            SINGDAT(:,1) = DSS_rmCommCond(SINGDAT(:,1),[],[2 1 3],rmCommCond,'off');
            if ~isempty(SINGDAT{1,2}); SINGDAT(:,2) = DSS_rmCommCond(SINGDAT(:,2),[],[2 1 3],rmCommCond,'off'); end
            for i = 1:no_cond
                DATs{s,1}{s2,1}{i,1} = mean(SINGDAT{i,1},3);
                if ~isempty(SINGDAT{i,2}); DATs{s,1}{s2,1}{i,2} = mean(SINGDAT{i,2},3);end
            end
        end
        
        if singtrial; save([sublist{s,s2} 'singdat_' sublist_side{s} '.mat'],'SINGDAT'); end
    end
end

if flip 
    suff = '_flip';
else
    suff='';
end
save([select '-' TFmethod savenametype typehand suff '_data.mat'],'-v7.3');

