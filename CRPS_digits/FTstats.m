function stat = FTstats(statmode,subjinfo,condlist,condcont,latency,frequency,cov,filepath,varargin)

param = finputcheck(varargin, {
    'alpha' , 'real' , [], 0.05; ...
    'numrand', 'integer', [], 1000; ...
    'ttesttail', 'integer', [-1 0 1], 0; ...
    'test_gfporchan', 'string', {}, 'off';...
    'singlesource', 'string', {'on' 'off'}, 'off';...
    'testmean', 'string', {'on' 'off'}, 'off';...
    'testlat', 'string', {'on' 'off'}, 'off';...
    'plotdata', 'cell', [], cell(1);...
    'timeshift' , 'real' , [], 0; ...
    'peakdef','cell' , [], cell(1); ...
    'gfpbasecorrect','real' , [], 0; ...
    'freqres','real' , [], 2; ...
    'numERPcomp','real' , [], []; ...
    'eventtypes','cell' , [], cell(1); ...
    'use_etype','real' , [], 0; ...
    });

timeshift = param.timeshift; %seconds

if iscell(frequency)
    freqsall = [];
    for f = 1:length(frequency)
        freqsall = [freqsall frequency{f}(1):param.freqres:frequency{f}(2)];
    end
    frequency = unique(freqsall);
end

if strcmp(statmode,'trial') && ischar(subjinfo)
    %%%% perform single-trial statistics
    subjlist = {subjinfo};
    subjcond = condlist;
    
elseif strcmp(statmode,'trial') && iscell(subjinfo)
    %%%% perform single-trial statistics
    subjlist = subjinfo;
    subjcond = condlist;
    
elseif strcmp(statmode,'cond') %&& length(subjinfo) == 2
    %%%% perform within-subject statistics
    try 
        if isstruct(subjinfo{1,1}{1,1}{1,1}{1,1})
            subjlist_temp = vertcat(subjinfo{:}); % creates a single column where Nrows=Nsubjects
            subjlist_temp = horzcat(subjlist_temp{:})'; % creates one column per condition
            diff = length(subjlist_temp{1,1});
            if diff==2
                subjlist = horzcat(subjlist_temp{:,1})';
                subjlist(:,3:4) = horzcat(subjlist_temp{:,2})';
            else
                subjlist = horzcat(subjlist_temp{:,1})';
                subjlist(:,2) = horzcat(subjlist_temp{:,2})';
            end
        end
    catch 
        subjlist1 = subjinfo{1};
        subjlist2 = subjinfo{2};
        subjlist = cat(2,subjlist1,subjlist2);
    end
    subjcond = repmat(condlist,length(subjlist),1);
    
    
elseif strcmp(statmode,'condtrial') && length(subjinfo) == 1
    %%%% perform within-subject statistics
    subjlist = subjinfo{1};
    subjcond = repmat(condlist,length(subjlist),1);
    
elseif strcmp(statmode,'subj') && iscell(subjinfo) && length(subjinfo) == 2 || length(subjinfo) == 4
    ucond = unique(condlist);
    subjlist=[];
    subjcond = [];
    numsubjgrp = [];
    for uc = 1:length(ucond)
        repcondi = find(strcmp(condlist,ucond(uc)));
        subjRC=[];
        condRC=[];
        for rc = 1:length(repcondi)
            subjRC = cat(2,subjRC,subjinfo{repcondi(rc)});
            condRC = cat(2,condRC, repmat(condlist(repcondi(rc)),length(subjinfo{repcondi(rc)}),1));
        end
        subjlist = cat(1,subjlist,subjRC);
        subjcond = cat(1,subjcond,condRC);
        numsubjgrp = [numsubjgrp size(subjRC,1)];
    end
    condlist = ucond;
    
elseif strcmp(statmode,'subj_corr') && iscell(subjinfo) && ~isempty(cov)
    if iscell(subjinfo{1,1})
        subjlist = subjinfo{:};
    else
        subjlist = subjinfo;
    end
    subjlist(:,2) = repmat({'cov'},size(subjlist,1),1);
    subjcond = condlist;
elseif strcmp(statmode,'corr') && iscell(subjinfo) && ~isempty(cov)
    if iscell(subjinfo{1,1})
        subjlist = subjinfo{:};
    else
        subjlist = subjinfo;
    end
    subjcond = condlist;
else
    error('Invalid combination of statmode and subjlist!');
end

numsubj = size(subjlist,1);
numdat = size(subjlist,2);

conddata = cell(numsubj,numdat);
tldata = cell(numsubj,numdat);

%% load and prepare individual subject datasets
for s = 1:numsubj
    %     reref[-200 0];
%     bcwin = erence
%     EEG = rereference(EEG,1);
    
%     %%%% baseline correction relative to 5th tone
%     bcwin = bcwin+(timeshift*1000);
%     EEG = pop_rmbase(EEG,bcwin);
%     %%%%
    
    for c = 1:numdat
        isEEG=0;
        if size(subjlist,2)>1
            if isstruct(subjlist{s,c})
                EEG = subjlist{s,c};
                isEEG=1;
            elseif strfind(subjlist{s,c},'spm')
                Ds = spm_eeg_load(fullfile(pwd,subjlist{s,c}));
                conddata{s,c} = spm2fieldtrip(Ds);
            elseif strcmp(subjlist{s,c},'cov')==1
                conddata{s,c} =conddata{s,1};
                if strfind(subjlist{s,1},'spm')
                    for t=1:size(conddata{s,c}.trial,2)
                        conddata{s,c}.trial{1,t} = cov(s)*ones(size(conddata{s,c}.trial{1,t}));
                    end
                else
                    conddata{s,c}.data(:) = cov(s);
                end
            else
                EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,c}), 'filepath', filepath);
                isEEG=1;
            end
        else 
            EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,1}), 'filepath', filepath);
            isEEG=1;
        end
        
        if isEEG
            % THIS ASSUMES THAT ALL DATASETS HAVE SAME NUMBER OF ELECTRODES
            if s == 1 && c == 1
                chanlocs = EEG.chanlocs;
            end

            %---split EEG into conditions based on trial indices ---%
            if strcmp(statmode,'condtrial')
                condevents = EEG2condevents(EEG,'STIM', param.eventtypes,param.use_etype); % currently a tailored function to CRPS digit perception study
                %no_cond = length(unique(condevents));
                EEG = pop_selectevent(EEG, 'event', find(condevents==str2double(condlist{c})));
            end

            conddata{s,c} = EEG;

            if isempty(frequency)
                % de-noise ERP
                if ~isempty(param.numERPcomp)
                    conddata{s,c}.data = DSSonERP(conddata{s,c}.data,[2 1 3],[],param.numERPcomp,[],'off'); % maximise repeatability
                end
            end
        end

%         if (strcmp(statmode,'trial') || strcmp(statmode,'cond')) && c == numcond
%             if conddata{s,1}.trials > conddata{s,2}.trials
%                 fprintf('Equalising trials in condition %s.\n',subjcond{s,1});
%                 randtrials = randperm(conddata{s,1}.trials);
%                 conddata{s,1} = pop_select(conddata{s,1},'trial',randtrials(1:conddata{s,2}.trials));
%             elseif conddata{s,2}.trials > conddata{s,1}.trials
%                 fprintf('Equalising trials in condition %s.\n',subjcond{s,2});
%                 randtrials = randperm(conddata{s,2}.trials);
%                 conddata{s,2} = pop_select(conddata{s,2},'trial',randtrials(1:conddata{s,1}.trials));
%             end
%         end
    end
end

%% prepare for fieldtrip statistical analysis
cfg = [];
cfg.keeptrials = 'yes';
cfg.feedback = 'textbar';

for s = 1:size(conddata,1)
    for c = 1:size(conddata,2)
        if ~strcmp(param.test_gfporchan,'off') && ~strcmp(param.test_gfporchan,'gfp') && (strcmp(statmode, 'condtrial') || strcmp(statmode, 'cond') || strcmp(statmode,'subj') || strcmp(statmode,'corr')) && ~any(strcmp('cfg',fieldnames(conddata{s,c})))
            if isempty(frequency)
                tldata{s,c} = ft_timelockanalysisCAB(cfg, convertoft(selectchan(conddata{s,c},param.test_gfporchan,'mean')));
            else
                tldata{s,c}= FTfreqanalysis(convertoft(selectchan(conddata{s,c},param.test_gfporchan)),frequency,conddata{s,c}.times/1000,conddata{s,c}.times(conddata{s,c}.times<0)/1000,1);
            end
        elseif strcmp(param.test_gfporchan,'gfp') && (strcmp(statmode, 'condtrial') || strcmp(statmode, 'cond') || strcmp(statmode,'subj') || strcmp(statmode,'corr')) && ~any(strcmp('cfg',fieldnames(conddata{s,c})))
            if isempty(frequency)
                tldata{s,c} = ft_timelockanalysisCAB(cfg, convertoft(convertogfp(conddata{s,c},param.gfpbasecorrect)));
            else
                tldata{s,c}= FTfreqanalysis(convertoft(conddata{s,c}),frequency,conddata{s,c}.times/1000,conddata{s,c}.times(conddata{s,c}.times<0)/1000,1);
            end
        elseif strcmp(param.test_gfporchan,'gfp') && strcmp(statmode,'subj_corr') && ~any(strcmp('cfg',fieldnames(conddata{s,c})))
            if c==1
                if isempty(frequency)
                    tldata{s,c} = ft_timelockanalysisCAB(cfg, convertoft(convertogfp(conddata{s,c},param.gfpbasecorrect)));
                else
                    tldata{s,c}= FTfreqanalysis(convertoft(conddata{s,c}),frequency,conddata{s,c}.times/1000,conddata{s,c}.times(conddata{s,c}.times<0)/1000,1);
                end
            elseif c==2
                if isempty(frequency)
                    tldata{s,c} = ft_timelockanalysisCAB(cfg, convertoft(meanchan(conddata{s,c})));
                else
                    error('subj_corr not set up for freq analysis');
                end
            end
        else
            if any(strcmp('cfg',fieldnames(conddata{s,c})))
                if isempty(frequency)
                    tldata{s,c} = ft_timelockanalysisCAB(cfg, conddata{s,c});
                else
                    tldata{s,c}= FTfreqanalysis(conddata{s,c},frequency,conddata{s,c}.times/1000,conddata{s,c}.times(conddata{s,c}.times<0)/1000,1);
                end
            else
                if isempty(frequency)
                    tldata{s,c} = ft_timelockanalysisCAB(cfg, convertoft(conddata{s,c}));
                else
                    tldata{s,c}= FTfreqanalysis(convertoft(conddata{s,c}),frequency,conddata{s,c}.times/1000,conddata{s,c}.times(conddata{s,c}.times<0)/1000,1);
                end
            end
        end
        if isempty(frequency)
            
            % if multiple latencies and peak definitions are specified for
            % ERP
            if iscell(latency) && iscell(param.peakdef)
                lats = [latency{:}];
                lats = sort(unique([lats{:}]));
                tldata{s,c}.trial = tldata{s,c}.trial(:,:,lats);
                tldata{s,c}.time = tldata{s,c}.time(lats);
                tldata{s,c}.avg = squeeze(mean(tldata{s,c}.trial,1));
            end
        else
            if strcmp(param.test_gfporchan,'gfp')
                tldata{s,c}=freqconvertogfp(tldata{s,c},param.gfpbasecorrect);
            end
            % if multiple latencies and peak definitions are specified
            if iscell(latency) && iscell(param.peakdef)
                lats = [latency{:}];
                lats = sort(unique([lats{:}])); % must use all lat from all freq
                if strcmp(param.test_gfporchan,'gfp')
                    tldata{s,c}.powspctrm = tldata{s,c}.powspctrm(:,:,lats);
                else 
                    tldata{s,c}.powspctrm = tldata{s,c}.powspctrm(:,:,:,lats);
                end 
                tldata{s,c}.time = tldata{s,c}.time(lats);
            end
            if length(frequency)>1 && length(frequency)~=length(tldata{s,c}.freq)
                fi = dsearchn(tldata{s,c}.freq',frequency')';
                tldata{s,c}.freq = tldata{s,c}.freq(fi);
                if strcmp(param.test_gfporchan,'gfp')
                    tldata{s,c}.powspctrm = tldata{s,c}.powspctrm(fi,:);
                else
                    tldata{s,c}.powspctrm = tldata{s,c}.powspctrm(:,:,fi,:);
                end
            end
        end
    end
end

elec = tldata{s,c}.elec;
save FT_layout elec

%% perform fieldtrip statistics
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.correctm = 'cluster';
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the permutation distribution.

cfg.tail = param.ttesttail;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = param.ttesttail;
if param.ttesttail == 0
    cfg.alpha = param.alpha/2;               % alpha level of the permutation test
else
    cfg.alpha = param.alpha;
end
cfg.clusteralpha = param.alpha;         % alpha level of the sample-specific test statistic that will be used for thresholding

cfg.numrandomization = param.numrand;      % number of draws from the permutation distribution

if strcmp(param.test_gfporchan,'off') && strcmp(param.singlesource,'off')
    % prepare_neighbours determines what sensors may form clusters
    cfg_neighb.method    = 'distance';
    cfg_neighb.neighbourdist = 4;
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb,convertoft(conddata{1,1}));
else
    cfg.neighbours = [];
    cfg.minnbchan = 0;               % minimum number of neighborhood channels that is required for a selected
end

if strcmp(statmode,'trial')
    
    if isempty(cov)
        error('covariate required for single-trial analysis');
    cfg.minnbchan = 2;               % minimum number of neighborhood channels that is required for a selected

    end
    
    if length(unique(cov))==2
        testtype = 'indepsamplesT';
    elseif length(unique(cov))>2
        testtype = 'correlationT';
    end
    
    cfg_ga = [];
    %cfg_ga.keepindividual = 'yes';
    %cond1data = ft_timelockgrandaverage(cfg_ga, tldata{:,1});
    %cond1data.avg = squeeze(mean(cond1data.individual,1))'; 
    if isempty(frequency)
        cond1data = ft_freqdescriptives(cfg_ga, tldata);
    else
        cond1data = tldata;
    end
    design(1,:) = cov;
    cfg.ivar     = 1;
    cfg.type = 'Spearman';
    
elseif strcmp(statmode,'cond') || strcmp(statmode, 'condtrial')
    
    %group statistics: we will perform within-subject comparison of subject
    %averages
    testtype = 'depsamplesT';
    
    if isempty(frequency)
        cfg_ga = [];
        cfg_ga.keepindividual = 'yes';
        
        if size(tldata,2)==4 && any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.avg = (tldata{ti,1}.avg - tldata{ti,2}.avg);
                tldata{ti,2}.avg = (tldata{ti,3}.avg - tldata{ti,4}.avg);
                tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial - tldata{ti,2}.trial);
                tldata{ti,2}.trial = NaN;%(tldata{ti,1}.trial - tldata{ti,2}.trial);
            end
            tldata(:,3:4) = [];
        end
        
        cond1data = ft_timelockgrandaverage(cfg_ga, tldata{:,1});
        cond2data = ft_timelockgrandaverage(cfg_ga, tldata{:,2});
        cond1data.avg = squeeze(mean(cond1data.individual,1));
        cond2data.avg = squeeze(mean(cond2data.individual,1));
    else
        cfg_ga = [];
        for i = 1:size(tldata,1)
            tldata{i,1} = ft_freqdescriptives(cfg_ga, tldata{i,1});
            tldata{i,2} = ft_freqdescriptives(cfg_ga, tldata{i,2});
            if size(tldata,2)==4 && any(condcont==-1)
                tldata{i,3} = ft_freqdescriptives(cfg_ga, tldata{i,3});
                tldata{i,4} = ft_freqdescriptives(cfg_ga, tldata{i,4});
            end
        end
        
        if size(tldata,2)==4 && any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.powspctrm = (tldata{ti,1}.powspctrm - tldata{ti,2}.powspctrm);
                tldata{ti,2}.powspctrm = (tldata{ti,3}.powspctrm - tldata{ti,4}.powspctrm);
                %tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial - tldata{ti,2}.trial);
            end
            tldata(:,3:4) = [];
        end
        
        cfg_ga.keepindividual = 'yes';
        cond1data = ft_freqgrandaverage(cfg_ga, tldata{:,1});
        cond2data = ft_freqgrandaverage(cfg_ga, tldata{:,2});
        cond1data.avg = squeeze(mean(cond1data.powspctrm,1));
        cond2data.avg = squeeze(mean(cond2data.powspctrm,1));
    end
    
    design = zeros(2,2*numsubj);
    design(1,:) = [ones(1,numsubj) 2*ones(1,numsubj)];
    design(2,:) = [1:numsubj 1:numsubj];
    cfg.ivar     = 1;
    cfg.uvar     = 2;
    
elseif strcmp(statmode,'subj')
    
    %group statistics: we will perform across-subject comparison of subject
    %averages
    testtype = 'indepsamplesT';
    
    cfg_ga = [];
    
    if isempty(frequency)
        cfg_ga.keepindividual = 'yes';

        if size(tldata,2) > 1 && ~any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.avg = (tldata{ti,1}.avg + tldata{ti,2}.avg) / 2;
                tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial + tldata{ti,2}.trial) / 2;
            end
        elseif size(tldata,2) > 1 && any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.avg = (tldata{ti,1}.avg - tldata{ti,2}.avg);
                tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial - tldata{ti,2}.trial);
            end
        end
        tldata(:,2) = [];

        cond1data = ft_timelockgrandaverage(cfg_ga, tldata{1:numsubjgrp(1),1});
        cond2data = ft_timelockgrandaverage(cfg_ga, tldata{numsubjgrp(1)+1:end,1});
        cond1data.avg = squeeze(mean(cond1data.individual,1))'; 
        cond2data.avg = squeeze(mean(cond2data.individual,1))';
    else
        cfg_ga = [];
        cfg_ga.keepindividual = 'yes';
        for i = 1:size(tldata,1)
            tldata{i,1} = ft_freqdescriptives(cfg_ga, tldata{i,1});
            if size(tldata,2) > 1
                tldata{i,2} = ft_freqdescriptives(cfg_ga, tldata{i,2});
            end
        end
        if size(tldata,2) > 1 && ~any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.powspctrm = (tldata{ti,1}.powspctrm + tldata{ti,2}.powspctrm) / 2;
                %tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial + tldata{ti,2}.trial) / 2;
            end
        elseif size(tldata,2) > 1 && any(condcont==-1)
            for ti = 1:size(tldata,1)
                tldata{ti,1}.powspctrm = (tldata{ti,1}.powspctrm - tldata{ti,2}.powspctrm);
                %tldata{ti,1}.trial = NaN;%(tldata{ti,1}.trial - tldata{ti,2}.trial);
            end
        end
        tldata(:,2) = [];

        cond1data = ft_freqgrandaverage(cfg_ga, tldata{1:numsubjgrp(1),1});
        cond2data = ft_freqgrandaverage(cfg_ga, tldata{numsubjgrp(1)+1:end,1});
        cond1data.avg = squeeze(mean(cond1data.powspctrm,1));
        cond2data.avg = squeeze(mean(cond2data.powspctrm,1));
    end
    
    design = zeros(1,numsubj);
    design(1,1:numsubjgrp(1)) = 1;
    design(1,numsubjgrp(1)+1:end)= 2;
    cfg.ivar  = 1;                   % number or list with indices, independent variable(s)
elseif strcmp(statmode,'subj_corr')
    
    %we will perform across-subject correlation with cov
    
    testtype = 'intersubcorr';
    
    if isempty(frequency)
        cfg_ga = [];
        cfg_ga.keepindividual = 'yes';
        cond1data = ft_timelockgrandaverage(cfg_ga, tldata{:,1});
        cond2data = ft_timelockgrandaverage(cfg_ga, tldata{:,2});
        
        zdata = zscore(reshape(cond1data.individual,size(cond1data.individual,1),size(cond1data.individual,2)*size(cond1data.individual,3)));
        cond1data.individual = reshape(zdata,size(cond1data.individual,1),size(cond1data.individual,2),size(cond1data.individual,3));
        zdata = zscore(reshape(cond2data.individual,size(cond2data.individual,1),size(cond2data.individual,2)*size(cond2data.individual,3)));
        cond2data.individual = reshape(zdata,size(cond2data.individual,1),size(cond2data.individual,2),size(cond2data.individual,3));
        
        cond1data.avg = squeeze(mean(cond1data.individual,1))'; 
        cond2data.avg = squeeze(mean(cond2data.individual,1))';
    else
        cfg_ga = [];
        for i = 1:size(tldata,1)
            tldata{i,1} = ft_freqdescriptives(cfg_ga, tldata{i,1});
            tldata{i,2} = ft_freqdescriptives(cfg_ga, tldata{i,2});
        end
        cfg_ga.keepindividual = 'yes';
        cond1data = ft_freqgrandaverage(cfg_ga, tldata{:,1});
        cond2data = ft_freqgrandaverage(cfg_ga, tldata{:,2});
        cond1data.avg = squeeze(mean(cond1data.powspctrm,1));
        cond2data.avg = squeeze(mean(cond2data.powspctrm,1));
    end
    
    design = zeros(2,2*numsubj);
    design(1,:) = [ones(1,numsubj) 2*ones(1,numsubj)];
    design(2,:) = [1:numsubj 1:numsubj];
    cfg.ivar     = 1;
    cfg.uvar     = 2;
    cfg.type = 'Spearman';
elseif strcmp(statmode,'corr')
    
    %we will perform across-subject correlation with cov
    
    testtype = 'correlationT';
    
    if isempty(frequency)
        cfg_ga = [];
        cfg_ga.keepindividual = 'yes';
        cond1data = ft_timelockgrandaverage(cfg_ga, tldata{:,1});
        cond1data.avg = squeeze(mean(cond1data.individual,1))'; 
    else
        cfg_ga = [];
        cfg_ga.keepindividual = 'yes';
        for i = 1:size(tldata,1)
            tldata{i,1} = ft_freqdescriptives(cfg_ga, tldata{i,1});
        end
        cond1data = ft_freqgrandaverage(cfg_ga, tldata{:,1});
        cond1data.avg = squeeze(mean(cond1data.powspctrm,1));
    end
    
    %design = zeros(1,1*numsubj);
    %design(1,:) = ones(1,numsubj);
    design(1,:) = cov;
    cfg.ivar     = 1;
    %cfg.uvar     = 1;
    cfg.type = 'Spearman';
end

% test mean or latency of the peak with a narrow latency window
if (strcmp(statmode,'cond') || strcmp(statmode,'subj') || strcmp(statmode,'subj_corr') || strcmp(statmode,'corr')) && (strcmp(param.testmean,'on') || strcmp(param.testlat,'on')) && isempty(frequency)
    cond1datatemp = cond1data;
    cond1data.time = [];
    for lat = 1:length(unique(param.peakdef))
        upeak = unique(param.peakdef);
        peakind = find(param.peakdef==upeak(lat));
        if iscell(latency)
            latNC = cond1datatemp.time(ismember([latency{:}],[latency{peakind}]));
            cond1data.time = [cond1data.time lat];
        else
            latNC = latency;
            cond1data.time = 0;
        end
        timeidx = cond1datatemp.time >= latNC(1)+timeshift & cond1datatemp.time <= latNC(end)+timeshift;
        for ind = 1:size(cond1datatemp.individual,1)
            for chan = 1:length(cond1datatemp.label)
                if strcmp(param.testmean,'on')
                    if isempty(frequency) || ~isfield(cond1datatemp,'powspctrm')
                        summinfo = mean(squeeze(cond1datatemp.individual(ind,chan,timeidx)));
                    else
                        summinfo = mean(squeeze(cond1datatemp.individual(ind,chan,freqsrange==frequency,timeidx)));
                    end
                elseif strcmp(param.testlat,'on')
                    if isempty(frequency) || ~isfield(cond1datatemp,'powspctrm')
                        summinfo = calclat(cond1datatemp.time(timeidx),squeeze(cond1datatemp.individual(ind,chan,timeidx))',50);
                    else
                        summinfo = calclat(cond1datatemp.time(timeidx),squeeze(cond1datatemp.individual(ind,chan,freqsrange==frequency,timeidx))',50);
                    end
                end
                cond1data.individual(ind,chan,freqsrange==frequency,lat) = summinfo(1);
            end
        end
    end
    cond1data.individual = cond1datatemp.individual(:,:,:,1:lat);
    cond1data.avg = squeeze(mean(cond1data.individual,1))';
    
    if exist('cond2data','var')
        for lat = 1:length(unique(param.peakdef))
            cond2datatemp = cond2data;
            cond2data.time = [];
            if iscell(latency)
                latNC = cond2datatemp.time(ismember([latency{:}],[latency{peakind}]));
                cond2data.time = [cond2data.time lat];
            else
                latNC = latency;
                cond2data.time = 0;
            end
            timeidx = cond2datatemp.time >= latNC(1)+timeshift & cond2datatemp.time <= latNC(end)+timeshift;
            for ind = 1:size(cond2datatemp.individual,1)
                for chan = 1:length(cond2datatemp.label)
                    if strcmp(param.testmean,'on')
                        if isempty(frequency) || ~isfield(cond2datatemp,'powspctrm')
                            summinfo = mean(squeeze(cond2datatemp.individual(ind,chan,timeidx)));
                        else
                            summinfo = mean(squeeze(cond2datatemp.individual(ind,chan,freqsrange==frequency,timeidx)));
                        end
                    elseif strcmp(param.testlat,'on')
                        if isempty(frequency) || ~isfield(cond2datatemp,'powspctrm')
                            summinfo = calclat(cond2datatemp.time(timeidx),squeeze(cond2datatemp.individual(ind,chan,timeidx))',50);
                        else
                            summinfo = calclat(cond2datatemp.time(timeidx),squeeze(cond2datatemp.individual(ind,chan,freqsrange==frequency,timeidx))',50);
                        end
                    end
                    cond2data.individual(ind,chan,freqsrange==frequency,lat) = summinfo(1);
                end
            end
        end
        cond2data.individual = cond2data.individual(:,:,:,1:lat);
        cond2data.avg = squeeze(mean(cond2data.individual,1))';
    end
end

if iscell(latency)
    latwin = [cond1data.time(1) cond1data.time(end)];
else
    latwin = latency;
end

cfg.design = design;            
cfg.statistic = testtype;
cfg.latency = latwin + timeshift;            % time interval over which the experimental conditions must be compared (in seconds)
cfg.feedback = 'textbar';
fprintf('\nComparing conditions using %d-tailed %s test\nat alpha of %.2f between %.2f-%.2f sec.\n\n', param.ttesttail, testtype, param.alpha, latwin);

if isempty(param.plotdata{1,1})
    diffcond = cond1data;
    diffcond.cond1avg = cond1data.avg;
    if exist('cond2data','var')
        diffcond.cond2avg = cond2data.avg;
        diffcond.avg = cond1data.avg - cond2data.avg;
    end
elseif ~isempty(param.plotdata{1,1}) 
    diffcond = tlplotdata(param.plotdata);
    diffcond.cond1avg = diffcond.avg;
    if exist('cond2data','var')
        diffcond.cond2avg = diffcond.avg;
    end
end

if ~isempty(frequency)
    diffconddimord = 'chan_freq_time';
else
    diffconddimord = 'chan_time';
end

if isempty(frequency)
    if isfield(diffcond,'cond2avg')
        [stat] = ft_timelockstatistics(cfg, cond1data, cond2data);
    else
        [stat] = ft_timelockstatistics(cfg, cond1data);
    end
else
    %cfg.frequency = frequency;
    if isfield(diffcond,'cond2avg')
        [stat] = ft_freqstatistics(cfg, cond1data, cond2data);
    else
        [stat] = ft_freqstatistics(cfg, cond1data);
    end
end

if iscell(subjinfo) 
    if ~iscell(subjinfo) 
        subjname = subjinfo{1,1};
        subjname = subjname(1:length(subjname)-3);
    else
        subjname = 'grpanalysis';
    end
end
if isnumeric(subjinfo)
    subjname = num2str(subjinfo);
end
if length(condlist)==2
    save2file = sprintf('%s_%s_%s-%s.mat',statmode,subjname,condlist{1},condlist{2});
else
    save2file = sprintf('%s_%s_%s.mat',statmode,subjname,condlist{1});
end
if ~exist('chanlocs','var'); load chanlocs; end;
stat.chanlocs = chanlocs; 
stat.cfg = cfg;
stat.condlist = condlist;
stat.diffcond = diffcond;
stat.diffconddimord = diffconddimord;
stat.timeshift = timeshift;
stat.statmode = statmode;
stat.subjinfo = subjname;

%plotclusters(stat);
%plotcorr(stat,cond1data,cov);

%if nargout == 0
    save(save2file, 'stat');
%end

function EEG = convertogfp(EEG,gfpbasecorrect)
EEG.nbchan = 1;
EEG.trials = 1;
EEG.icachansind = 1;
EEG.chanlocs = EEG.chanlocs(1);
%[junk, gfp] = evalc('eeg_gfp(mean(EEG.data,3)'')'''); % Lehmann's original GFP
gfp = squeeze(std(mean(EEG.data,3),1));
if gfpbasecorrect; gfp = rmbase(gfp,[],1:find(EEG.times == 0));end;
EEG.data = gfp;

function EEG = freqconvertogfp(EEG,gfpbasecorrect)
EEG.label = EEG.label(1);
meanpow = squeeze(nanmean(EEG.powspctrm,1));
gfp = nanstd(reshape(meanpow,size(meanpow,1),size(meanpow,2)*size(meanpow,3)),1);
gfp = reshape(gfp,1,size(meanpow,2),size(meanpow,3));
if gfpbasecorrect; gfp = rmbase(gfp,[],1:find(EEG.times == 0));end;
EEG.powspctrm = gfp;
EEG.dimord = 'chan_freq_time';

function EEG = meanchan(EEG)
EEG.nbchan = 1;
EEG.trials = 1;
EEG.icachansind = 1;
EEG.chanlocs = EEG.chanlocs(1);
mc = squeeze(mean(mean(EEG.data,3),1));
EEG.data = mc;

function EEG = selectchan(EEG,chan,meanopt)
EEG.nbchan = 1;
EEG.icachansind = 1;
chani = strcmp({EEG.chanlocs.labels},chan);
EEG.chanlocs = EEG.chanlocs(chani);
EEG.data = EEG.data(chani,:,:);
if strcmp(meanopt,'mean')
    EEG.data = mean(EEG.data,3);
    EEG.trials = 1;
end

function estlat = calclat(times,data,pcarea)
%estlat = sum(abs(data));

totalarea = sum(abs(data));
pcarea = totalarea * (pcarea/100);

curarea = 0;
for t = 1:length(data)
    curarea = curarea + abs(data(t));
    if curarea >= pcarea
        estlat = times(t);
        return
    end
end
estlat = times(end);

function plotdata = tlplotdata(files)
subjlist = files;
subjcond = {'1'};
numsubj = size(subjlist,1);
numcond = size(subjcond,2);
conddata = cell(numsubj,numcond);
plotdata = cell(numsubj,numcond);
for s = 1:numsubj
    for c = 1:numcond
        if strfind(subjlist{s,c},'spm')
            D = spm_eeg_load(fullfile(pwd,subjlist{s,c}));
            conddata{s,c} = spm2fieldtrip(D);
            clear D classD
        else
            EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,c}), 'filepath', filepath);
            if s == 1 && c == 1
                chanlocs = EEG.chanlocs;
            end
            EEG.icachansind = 1:length(chanlocs);
            conddata{s,c} = EEG;
            clear EEG classEEG
        end
    end
end

% prepare for fieldtrip statistical analysis
cfg = [];
cfg.keeptrials = 'yes';
cfg.feedback = 'textbar';

for s = 1:size(conddata,1)
    for c = 1:size(conddata,2)
        if any(strcmp('cfg',fieldnames(conddata{s,c})))
            plotdata{s,c} = conddata{s,c};
        else
            plotdata{s,c} = convertoft(conddata{s,c});
        end
        conddata{s,c}=[];
        plotdata{s,c} = ft_timelockanalysisCAB(cfg, plotdata{s,c});
        plotdata{s,c}.trial = single(plotdata{s,c}.trial);
    end
end

cfg_ga = [];
cfg_ga.keepindividual = 'yes';
plotdata = ft_timelockgrandaverage(cfg_ga, plotdata{:});
plotdata.avg = squeeze(mean(plotdata.individual,1)); 

function EEG = rmerp(EEG)
trialdata = EEG.data;
trialdata2=trialdata;
erpdata = double(squeeze(mean(trialdata,3)));
for t = 1:size(trialdata,3)
    trialdata2(:,:,t) = trialdata(:,:,t)-erpdata; % subtract out ERP for induced activity
end
EEG.data = trialdata2;
