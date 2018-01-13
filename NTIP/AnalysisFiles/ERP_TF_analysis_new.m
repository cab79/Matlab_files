%% ERP, Frequency, Time-Frequency and Coherence analysis
% To be run after data can be cleaned and saved as ALLEEG structures (one
% per subject) - i.e after preprocess part 2
% Before running this, check data quality using "Plot_ALLEEG.m"
% Requires Fieldtrip to be in the Matlab path
% Save ALL SUBJECTS' data in one huge data file

clear all
close all
dbstop if error
%% PATHS AND DATA
% path to cleaned EEG data
filepath = 'C:\Data\NTIP\Preprocessed';
cd(filepath);
% generic file suffix
filesuff = ('_cleaned.set');
filepref ='NTIP_';
% chanlocs file path
load('C:\Data\NTIP\chanlocs.mat');

% load .xlsx file containing columns named 'Participant_ID', 'Group', and
% any covariates of interest
pdatfile = 'C:\Data\NTIP\Participant_data.xlsx';

%% SETTINGS

% Here select analysis by selecting one line at a time:
%select = 'TF'; TFmethod = '-FT'; % or, '-EL' % recommend FT (fieldtrip)
select = 'Freq'; TFmethod = '-FT'; % (fieldtrip)
%select = 'Coh'; TFmethod = '-FT'; % Coherence (fieldtrip)
%select = 'ERP'; TFmethod = ''; % if ERP, leave TF blank

timebin = [-0.2 0.399]; % time window
basebin = [-0.2 0]; % baseline window

% select frequency range and resolution (if freq analysis)
freqsrange = [8:2:12]; % example of keeping freq resolution high (2Hz) for TF analysis
%freqsrange = [0 4; 4 8; 8 13; 13 30; 30 40]; % example of choosing canonical windows for Freq analysis

% Bootstraps repetitions (for coherence analysis only)
bootrep = 50;

% number of cycles (larger the better but may crash if too big)
ncycles=3;

% files to load
filetypes = {'10HC','10HO','10LC','10LO','1HC','1HO','1LC','1LO'};

% stimtypes to include in analysis
eventtypes = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8', 'S  9'}; % stim labels
use_etype = [1 2 3 4 5 6 7 8 9]; % index of labels to include
no_cond = length(use_etype);

% save single trial data? (large files!)
singtrial=0;

% Settings for DSS (recommend not using this)
DSS_ERP=0; 
numERPcomp=4;
rmCommCond=0;

% Apply CSD (leave as 0 unless you know what you are doing!)
applyCSD=0; 
% path to CSD montage (only needed if using CSD)
if strcmp(select,'ERP') && applyCSD; load M:\Matlab\Matlab_files\CRPS_digits\CSDmontage_92.mat; end



%% RUN

cd(filepath)

% load participant info and identify the relevant columns
[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),'Group'));
sub_col = find(strcmp(pdata(1,:),'Subject'));
inc_col = find(strcmp(pdata(1,:),'Include'));

% identify number of groups of participants and find indices of
% participants in each group
Ngrp = length(unique([pdata{2:end,grp_col}]));
SubInd = cell(Ngrp,1);
Subs = [];
gn=0;
for g = unique([pdata{2:end,grp_col}])
    gn=gn+1;
    inc_idx = find(cellfun(@(x) x==1, pdata(2:end,inc_col), 'UniformOutput', 1));
    grp_idx = find(cellfun(@(x) x==g, pdata(2:end,grp_col), 'UniformOutput', 1));
    SubInd{gn,1} = intersect(inc_idx,grp_idx);
    Nsub(gn,1) = length(SubInd{gn,1});
end

% make some corrections to the subject IDs if needed to make life easier
% later
subjlists={};
for g = 1:Ngrp
    subgrp={};
    s2=0;
    for s = 1:Nsub(g)
        %if isnan(group(s,g))
        %    continue
        %end
        s2 = s2+1;
        subj = pdata{1+SubInd{g}(s),sub_col};
        if isnumeric(subj)
            subj = num2str(subj);
        end
        subgrp{s2,1} = subj;
    end
    subjlists{g,1} = subgrp;
end
grplist = 1:Ngrp; % index of groups
el = length(chanlocs); % number of channels

% ignore these bits
% for each group, was left or right arm stimulated? Leave as '' if
% irrelevant
sublist_side = {''}; 
% for each group, how should they be labelled?
sublist_grp = {'x'}; 

savenametype='';
for i = use_etype
    savenametype = [savenametype '_' eventtypes{1,i}(end)];
end

%el = 62;

% create DAT structures for analysed data
subjects = subjlists(grplist);
%DATs = subjects;

if strcmp(select,'ERP')
%    dp = 250;
    nf = 1;
%    times = -0.2:0.004:0.796;
elseif strcmp(select,'TF')
%    dp = 181;
    nf = length(freqsrange);
%    intimes = -0.2:0.004:0.796;
end

% create DATs and grand average
sublist = {};
for g = 1:length(subjects)
    for s = 1:length(subjects{g,1}) 
        subj = subjects{g,1}{s,1};
        subjpath = fullfile(filepath);
        for f = 1:length(filetypes)
            fname = [filepref subj '_' filetypes{f} filesuff];
            EEG = pop_loadset(fname,subjpath);
            
            % Remove ECG
            ecgchan = strcmp({EEG.chanlocs.labels},'ECG');
            if any(ecgchan)
                EEG = pop_select(EEG,'nochannel',find(ecgchan));
            end
            
            dsize = size(EEG.data);
            intimes=EEG.times;

            % data cleaning and ERP - not recommended
            %EEGrep=EEG; % most repeatable over trials
            %EEGrep.data = DSSonERP(EEGrep.data,[2 1 3],condevents,[4],[],'off'); % maximise repeatability
            %pause
            %close all

            DATs{g,1}{s,1}{f,1} = cell(no_cond,1);
            SINGDAT = cell(1);
            EEGall=EEG;
            for st = 1:no_cond
                EEG = pop_selectevent(EEGall,'type',eventtypes{use_etype(st)});
                EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
                eegdata = [];
                eegdatasing = [];
                inddata = [];
                inddatasing = [];
                evdata = [];
                evdatasing = [];
                itcdata = [];
                matrix = [];
                bootmat=[];
                if length(EEG.epoch)>1
                   evdatasing = EEG.data;
                   times=intimes;
                    if strcmp(select,'ERP') && applyCSD; 
                        [nchan,nsamp,ntrial] = size(evdatasing);
                        evdatasing = reshape(evdatasing,nchan,nsamp*ntrial);
                        evdatasing=CSD(evdatasing,G,H); 
                        evdatasing = reshape(evdatasing,nchan,nsamp,ntrial);
                    end
                    if DSS_ERP; evdatasing = DSSonERP(evdatasing,[2 1 3],[],[numERPcomp],[],'off');end; % maximise repeatability 
                    evdata = double(mean(evdatasing,3));
                    if strcmp(select,'TF') || strcmp(select,'Freq') || strcmp(select,'Coh')
                        trialdata= EEG.data;
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
                                [eegdata(:,:,e),eegdatasing(:,:,:,e),itcdata(:,:,e),~,times,freqs,~,~] = newtimef_singletrial(datvec,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                                [eegdata2(:,:,e),eegdatasing2(:,:,:,e),~,~,times,freqs,~,~] = newtimef_singletrial(datvec2,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'freqs',[freqsrange(1) freqsrange(end)],'nfreqs',nf,'baseline',0,'scale','log','plotersp','off','plotitc','off','plotphase','off');
                                %[eegdata(:,:,e),eegdatasing(:,:,:,e),~,~,times,freqs] = timef_singletrial(datvec,size(EEG.data,2),[-200 800],EEG.srate,[1 0.5],'detrend','off','plotersp','off','plotitc','off','plotphase','off');  
                            end
                            inddata = permute(eegdata2,[3 2 1]);
                            evdata = permute(eegdata-eegdata2,[3 2 1]);
                            itcdata = permute(itcdata,[3 2 1]);
                            inddatasing = permute(eegdatasing2,[4 2 3 1]);
                            evdatasing = permute(eegdatasing-eegdatasing2,[4 2 3 1]);
                        elseif strcmp(TFmethod,'-FT')
                            FTEEG = EEG;%pop_selectevent(EEG,'event',condidx);
                            FTEEG2 = FTEEG;
                            FTEEG2.data = trialdata2;
                            FTEEG2=convertoft(FTEEG2);
                            FTEEG=convertoft(FTEEG);

                            % freq analysis
                            FTEEG2 = FTfreqanalysis(FTEEG2,select,freqsrange,intimes,basebin,ncycles);
                            inddatasing = permute(FTEEG2.powspctrm,[2 4 1 3]);
                            inddata = squeeze(nanmean(inddatasing,3));
                            act_freq = FTEEG2.freq;
                            disp(['actual frequencies: ' num2str(act_freq)])

                            if strcmp(select,'Coh')
                                % coherence matrices
                                matrix=zeros(size(freqsrange,1),EEG.nbchan,EEG.nbchan); 
                                coh = zeros(EEG.nbchan,EEG.nbchan);
                                cohboot = zeros(EEG.nbchan,EEG.nbchan);

                                wpli = ft_connectivity_wpli(FTEEG2.crsspctrm,'debias',true,'dojack',false);
                                for f = 1:size(freqsrange,1)
                                    fprintf('f %d',f);
                                    [M, bstart] = min(abs(FTEEG2.freq-freqlist(f,1)));
                                    [M, bend] = min(abs(FTEEG2.freq-freqlist(f,2)));
                                    coh(:) = 0;
                                    coh(logical(tril(ones(size(coh)),-1))) = max(wpli(:,bstart:bend),[],2);
                                    coh = tril(coh,1)+tril(coh,1)';
                                    matrix(f,:,:) = coh;
                                end

                                for nboot = 1:bootrep
                                    fprintf('nboot %d',nboot);
                                    for tr = 1:length(FTEEG2.trial)
                                        for ele = 1:size(FTEEG2.trial{tr},1);
                                            trial = FTEEG2.trial{tr};
                                            surrEEG=[];
                                            surrEEG = [phaseran(squeeze(trial(tr,:)'),1); 0]';
                                            %surrEEG = [phaseran(EEG.trial{1}',1); zeros(1, length(EEG.label))]';
                                            trial(tr,:) = surrEEG;
                                        end
                                        FTEEG2.trial{tr} = trial;
                                    end

                                    FTEEG2 = FTfreqanalysis(FTEEG2,select,freqsrange,intimes,basebin,ncycles);
                                    wpli_boot = ft_connectivity_wpli(FTEEG2.crsspctrm,'debias',true,'dojack',false);

                                    for f = 1:size(freqlist,1)
                                        [M, bstart] = min(abs(FTEEG2.freq-freqlist(f,1)));
                                        [M, bend] = min(abs(FTEEG2.freq-freqlist(f,2)));

                                        cohboot(:) = 0;
                                        cohboot(logical(tril(ones(size(cohboot)),-1))) = max(wpli_boot(:,bstart:bend),[],2);
                                        cohboot = tril(cohboot,1)+tril(cohboot,1)';
                                        bootmat(f,:,:,nboot) = cohboot;
                                    end
                                    clear wpli_boot
                                end

                            elseif strcmp(select,'TF')
                                FTEEG = FTfreqanalysis(FTEEG,select,freqsrange,intimes,basebin,ncycles);
                                evdatasing = permute(FTEEG.powspctrm,[2 4 1 3])-inddatasing;
                                evdata = squeeze(nanmean(evdatasing,3));
                            end

                            times = FTEEG.time;
                            freqs = freqsrange;

                            if strcmp(select,'TF')
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
                end

                try
                    SINGDAT{st,1} = evdatasing;
                    SINGDAT{st,2} = inddatasing;
                    SINGDAT{st,3} = length(EEG.epoch);
                end
                
                if strcmp(select,'ERP')
                    DATs{g,1}{s,1}{f,1}{st,1} = evdata;
                end
                DATs{g,1}{s,1}{f,1}{st,2} = inddata;
                DATs{g,1}{s,1}{f,1}{st,3} = itcdata;
                DATs{g,1}{s,1}{f,1}{st,4} = matrix; % wpli coherence
                DATs{g,1}{s,1}{f,1}{st,5} = bootmat;
                DATs{g,1}{s,1}{f,1}{st,6} = length(EEG.epoch);
                DATs{g,1}{s,1}{f,1}{st,7} = act_freq;
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

            if singtrial; save([subj '_singdat.mat'],'SINGDAT'); end
        end
    end
end

save([select TFmethod savenametype '_data.mat'],'-v7.3');

%% average over conditions for each subject and file; then plot a frequency
load([select TFmethod savenametype '_data.mat']);
% settings
freqselect = 10;
condselect = 1;%[1:9];
exclchan = [29 30];
%block_reorder = 1:8; % no reordering
%filetype_subtract = [1 5; 2 6; 3 7; 4 8]; % subtract 10Hz v 1 Hz
%filetype_subtract = [1 2; 3 4; 5 6; 7 8]; % subtract Oddball v control
%filetype_subtract = [1 3; 2 4; 5 7; 6 8]; % subtract High vs Low freq

block_reorder = [4 1 8 5 7 6 3 4]; % re-order in time
filetype_subtract = [1 5; 2 6; 3 7; 4 8]; % subtract 10Hz v 1 Hz
participant_blockorder = [1 -1 1 -1]; % normal or reversed ordering

% run
close all
inclchan = 1:length(chanlocs);
inclchan(exclchan) = [];
meandat={};
f1=figure
f2=figure
gs=0;
gsf1=0;
gsf2=0;
dat_sub={};
for g = 1:length(subjects)
    for s = 1:length(subjects{g,1}) 
        gs=gs+1;
        if participant_blockorder(s)==-1
            block_order = fliplr(block_reorder);
        else
            block_order = block_reorder;
        end
        for f = 1:length(filetypes)
            dat = cat(3,DATs{g,1}{s,1}{block_order(f),1}{:,2});
            meandat{f,1} = mean(dat(:,:,condselect),3);
            act_freq=DATs{g,1}{s,1}{f,1}{1,7}
        end
        datmat = cat(3,meandat{:});
        maxdat = max(abs(datmat(:)));
        minmax = [-maxdat,maxdat];
        
        %% plot a frequency
        fn = dsearchn(act_freq',freqselect);
        for f = 1:length(meandat)
            gsf1=gsf1+1;
            dat = meandat{f}(:,fn);
            figure(f1)
            subplot(length(subjects{g,1}),length(meandat),gsf1)
            topoplot(dat(inclchan),chanlocs(inclchan),'maplimits',minmax);
            cbar('vert',0,minmax);
            title(['subject ' num2str(s) '  ' filetypes{f}]);
        end
        
        % calc subtraction
        for su = 1:length(filetype_subtract)
            dat1 = meandat{filetype_subtract(su,1)}(:,fn);
            dat2 = meandat{filetype_subtract(su,2)}(:,fn);
            dat_sub{su} = dat1-dat2;
        end
        
        %plot
        datmat = cat(3,dat_sub{:});
        maxdat = max(abs(datmat(:)));
        minmax = [-maxdat,maxdat];
        for su = 1:length(filetype_subtract)
            figure(f2)
            gsf2=gsf2+1;
            figure(f2)
            subplot(length(subjects{g,1}),length(filetype_subtract),gsf2)
            topoplot(dat_sub{su}(inclchan),chanlocs(inclchan),'maplimits',minmax);
            cbar('vert',0,minmax);
            title(['subject ' num2str(s) '  ' filetypes{filetype_subtract(su,1)} '-' filetypes{filetype_subtract(su,2)}]);
        end
    end
end


