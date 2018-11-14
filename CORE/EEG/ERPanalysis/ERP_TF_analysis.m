clear all
close all

filepath = 'C:\Data\CORE\eeg\ana\prep\cleaned\part4';
cd(filepath);
filesuff = ('_merged_cleaned.set');
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat');
%load('C:\Data\CORE\Preprocessed_100Hz\groupsublist.mat');
subjlists={};
sublist_grp = {'CORE'}; grplist = [1]; sublist_side = {'LR'}; expts = {'4'}; 
hand_nme = {'L','R'};

for g = 1:size(group,2)
    subgrp={};
    s2=0;
    for s = 1:size(group,1)
        if isnan(group(s,g))
            continue
        end
        s2 = s2+1;
        subj = num2str(group(s,g));
        if length(subj)<2
            subj = ['00' subj];
        elseif length(subj)<3
            subj = ['0' subj];
        end
        subgrp{s2,1} = [sublist_grp{g} subj];
    end
    subjlists{g,1} = subgrp;
end

el = length(chanlocs);

timebins = [-0.2 0.9; % for epoching, TSOT(2)
            -0.2 0.3]; % for epoching, TSOT(4)
basebins = [-0.2 0; % for epoching, TSOT(2)
            -0.05 0]; % for epoching, TSOT(4)
        
eventtypes = {'CNUM'; [1:24]};
use_etype = [1];

%eventtypes = cellfun(@num2str, num2cell(use_etype), 'UniformOutput', false);
%no_cond = 24;

%select = 'TF'; TFmethod = '-FT'; % or, '-EL'
select = 'ERP'; TFmethod = '';
singtrial=0;
numERPcomp=4;
DSS_ERP=0; 
rmCommCond=0;
applyCSD=0; 
if strcmp(select,'ERP') && applyCSD; load M:\Matlab\Matlab_files\CRPS_digits\CSDmontage_92.mat; end
%--------------%

savenametype='';
for i = use_etype
    savenametype = [savenametype '_' eventtypes{1,i}];
end

cd(filepath)

%el = 62;

subjects = subjlists(grplist);
DATs = subjects;

if strcmp(select,'ERP')
%    dp = 250;
    nf = 1;
%    times = -0.2:0.004:0.796;
elseif strcmp(select,'TF')
%    dp = 181;
    freqsrange = [2:2:40];
    nf = length(freqsrange);
%    intimes = -0.2:0.004:0.796;
end
    

for ex = 1:length(expts)
    
    if strcmp(expts{ex},'2')
        timebin = timebins(1,:);
        basebin = basebins(1,:);
    elseif strcmp(expts{ex},'4')
        timebin = timebins(2,:);
        basebin = basebins(2,:);
    end
    
    % create DATs and grand average
    sublist = {};
    for g = 1:length(subjects)
        for s = 1:length(subjects{g,1}) 
            subj = subjects{g,1}{s,1};
            EEG = pop_loadset([subj '_' expts{ex} filesuff],filepath);
            dsize = size(EEG.data);
            intimes=EEG.times;
            
            condevents = EEG2condevents(EEG,'STIM', eventtypes,use_etype);
            no_cond = length(unique(condevents));

            % data cleaning and ERP
            %EEGrep=EEG; % most repeatable over trials
            %EEGrep.data = DSSonERP(EEGrep.data,[2 1 3],condevents,[4],[],'off'); % maximise repeatability
            %pause
            %close all

            DATs{g,1}{s,1} = cell(no_cond,1);
            SINGDAT = cell(1);
            EEGall=EEG;
            for i = 1:no_cond
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
                   times=intimes;
                    if strcmp(select,'ERP') && applyCSD; 
                        [nchan,nsamp,ntrial] = size(evdatasing);
                        evdatasing = reshape(evdatasing,nchan,nsamp*ntrial);
                        evdatasing=CSD(evdatasing,G,H); 
                        evdatasing = reshape(evdatasing,nchan,nsamp,ntrial);
                    end
                    if DSS_ERP; evdatasing = DSSonERP(evdatasing,[2 1 3],[],[numERPcomp],[],'off');end; % maximise repeatability 
                    evdata = double(mean(evdatasing,3));
                    if strcmp(select,'TF')
                        trialdata= EEG.data(:,:,condidx);
                        trialdata2= trialdata;
                        for t = 1:size(trialdata,3)
                            trialdata2(:,:,t) = trialdata(:,:,t)-evdata; % subtract out ERP for induced activity
                        end
                        eegdata2 = [];
                        eegdatasing2 = [];
                        if strcmp(TFmethod,'-EL')
                            for e = 1:size(EEG.data,1)
                                datvec = reshape(squeeze(trialdata(e,:,:)),1,size(trialdata,2)*size(trialdata,3));
                                datvec2 = reshape(squeeze(trialdata2(e,:,:)),1,size(trialdata2,2)*size(trialdata2,3));
                                %[eegdata(:,:,e),x,y,times,freqs] = timef(datvec,size(EEG.data,2),[-200 800],EEG.srate,1,'detret','on','winsize',60,'plotersp','off','plotitc','off','plotphase','off');  
                                [eegdata(:,:,e),eegdatasing(:,:,:,e),itcdata(:,:,e),~,times,freqs,~,~] = newtimef_singletrial(datvec,size(EEG.data,2),timebin,EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                                [eegdata2(:,:,e),eegdatasing2(:,:,:,e),~,~,times,freqs,~,~] = newtimef_singletrial(datvec2,size(EEG.data,2),timebin,EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                                %[eegdata(:,:,e),eegdatasing(:,:,:,e),~,~,times,freqs] = timef_singletrial(datvec,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'detrend','off','plotersp','off','plotitc','off','plotphase','off');  
                            end
                            inddata = permute(eegdata2,[3 2 1]);
                            evdata = permute(eegdata-eegdata2,[3 2 1]);
                            itcdata = permute(itcdata,[3 2 1]);
                            inddatasing = permute(eegdatasing2,[4 2 3 1]);
                            evdatasing = permute(eegdatasing-eegdatasing2,[4 2 3 1]);
                        elseif strcmp(TFmethod,'-FT')
                            FTEEG = pop_select(EEG,'trial',condidx);
                            FTEEG2 = FTEEG;
                            FTEEG2.data = trialdata2;
                            FTEEG2=convertoft(FTEEG2);
                            FTEEG=convertoft(FTEEG);

                            % freq analysis
                            ncycles=3;
                            FTEEG = FTfreqanalysis(FTEEG,freqsrange,intimes,basebin,ncycles);
                            FTEEG2 = FTfreqanalysis(FTEEG2,freqsrange,intimes,basebin,ncycles);

                            inddatasing = permute(FTEEG2.powspctrm,[2 4 1 3]);
                            evdatasing = permute(FTEEG.powspctrm,[2 4 1 3])-inddatasing;
                            inddata = squeeze(nanmean(inddatasing,3));
                            evdata = squeeze(nanmean(evdatasing,3));
                            times = FTEEG.time;
                            freqs = freqsrange;

                            % remove baseline from average (baseline already
                            % removed from single trials)
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

                DATs{g,1}{s,1}{i,1} = evdata;
                DATs{g,1}{s,1}{i,2} = inddata;
                DATs{g,1}{s,1}{i,3} = itcdata;
                DATs{g,1}{s,1}{i,4} = length(condidx);
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
                    DATs{g,1}{s,1}{i,1} = mean(SINGDAT{i,1},3);
                    if ~isempty(SINGDAT{i,2}); DATs{g,1}{s,1}{i,2} = mean(SINGDAT{i,2},3);end
                end
            end

            if singtrial; save([subj '_singdat_' sublist_side{g} '_' expts{ex} '.mat'],'SINGDAT'); end
        end
    end

    save([select TFmethod savenametype '_' expts{ex} '_data.mat'],'-v7.3');
end

