function stat = comperp(statmode,subjinfo,condlist,latency,varargin)

loadpaths

param = finputcheck(varargin, {
    'alpha' , 'real' , [], 0.05; ...
    'numrand', 'integer', [], 1000; ...
    'ttesttail', 'integer', [-1 0 1], 0; ...
    'testgfp', 'string', {'on' 'off'}, 'off';...
    'testcnv', 'string', {'on' 'off'}, 'off';...
    'testlat', 'string', {'on' 'off'}, 'off';...
    });

timeshift = 0; %seconds

loadsubj

if strcmp(statmode,'trial') && ischar(subjinfo)
    %%%% perform single-trial statistics
    subjlist = {subjinfo};
    subjcond = condlist;
    
elseif strcmp(statmode,'cond') && isnumeric(subjinfo) && length(subjinfo) == 2
    %%%% perform within-subject statistics
    subjlist1 = subjlists{subjinfo(1)};
    subjlist2 = subjlists{subjinfo(2)};
    subjlist = cat(2,subjlist1,subjlist2);
    subjcond = repmat(condlist,length(subjlist),1);
    
elseif strcmp(statmode,'subj') && isnumeric(subjinfo) && length(subjinfo) == 2 || length(subjinfo) == 4
    %%%% perform across-subject statistics
    subjlist1 = subjlists{subjinfo(1)};
    subjlist2 = subjlists{subjinfo(2)};
    
    numsubj1 = length(subjlist1);
    numsubj2 = length(subjlist2);
    subjlist = cat(1,subjlist1,subjlist2);
    subjcond = cat(1,repmat(condlist(1),numsubj1,1),repmat(condlist(2),numsubj2,1));
    if length(condlist) == 4
        subjcond = cat(2,subjcond,cat(1,repmat(condlist(3),numsubj1,1),repmat(condlist(4),numsubj2,1)));
        subjlist3 = subjlists{subjinfo(3)};
        subjlist4 = subjlists{subjinfo(4)};
        subjlist = cat(2,subjlist,cat(1,subjlist3,subjlist4));
    end
else
    error('Invalid combination of statmode and subjlist!');
end

numsubj = length(subjlist);
numcond = size(subjcond,2);

conddata = cell(numsubj,numcond);
tldata = cell(numsubj,numcond);

%% load and prepare individual subject datasets
for s = 1:numsubj
    

    
    
    for c = 1:numcond
        EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,c}), 'filepath', filepath);
        EEG = sortchan(EEG);
        
        %     rereference
   
        EEG = pop_reref( EEG, [35 81] );
        
        %     %%%% baseline correction relative to 5th tone
    bcwin = [-200 0];
    bcwin = bcwin+(timeshift*1000);
    EEG = pop_rmbase(EEG,bcwin);


        % THIS ASSUMES THAT ALL DATASETS HAVE SAME NUMBER OF ELECTRODES
        if s == 1 && c == 1
            chanlocs = EEG.chanlocs;
        end
        
        conddata{s,c} = EEG;
        
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
        if strcmp(param.testgfp,'on') && (strcmp(statmode, 'cond') || strcmp(statmode,'subj'))
            tldata{s,c} = ft_timelockanalysis(cfg, convertoft(convertogfp(conddata{s,c})));
        else
            ftdata = convertoft(conddata{s,c});
            tldata{s,c} = ft_timelockanalysis(cfg, ftdata);
        end
    end
end

%% perform fieldtrip statistics
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo ('montecarlo')or ('analytic') Method to calculate the significance probability
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

if strcmp(param.testgfp,'off')
    % prepare_neighbours determines what sensors may form clusters
    cfg_neighb.method    = 'distance';
    cfg_neighb.neighbourdist = 4;
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb,convertoft(conddata{1,1}));
    cfg.minnbchan = 2;               % minimum number of neighborhood channels that is required for a selected

else
    cfg.neighbours = [];
    cfg.minnbchan = 0;               % minimum number of neighborhood channels that is required for a selected
end

if strcmp(statmode,'trial')
    
    %single-subject statistics: we will compare potentially different
    %number of trials in the two conditions for this subject.
    ttesttype = 'indepsamplesT';
    
    design = zeros(1,size(tldata{1}.trial,1) + size(tldata{2}.trial,1));
    design(1,1:size(tldata{1}.trial,1)) = 1;
    design(1,size(tldata{1}.trial,1)+1:end)= 2;
    cfg.ivar  = 1;                   % number or list with indices, independent variable(s)
    
    cond1data = tldata{1};
    cond2data = tldata{2};
    
elseif strcmp(statmode,'cond')
    
    %group statistics: we will perform within-subject comparison of subject
    %averages
    ttesttype = 'depsamplesT';
    
    cfg_ga = [];
    cfg_ga.keepindividual = 'yes';
    cond1data = ft_timelockgrandaverage(cfg_ga, tldata{:,1});
    cond2data = ft_timelockgrandaverage(cfg_ga, tldata{:,2});
    cond1data.avg = squeeze(mean(cond1data.individual,1));
    cond2data.avg = squeeze(mean(cond2data.individual,1));
    
    design = zeros(2,2*numsubj);
    design(1,:) = [ones(1,numsubj) ones(1,numsubj)+1];
    design(2,:) = [1:numsubj 1:numsubj];
    cfg.ivar     = 1;
    cfg.uvar     = 2;
    
elseif strcmp(statmode,'subj')
    
    %group statistics: we will perform across-subject comparison of subject
    %averages
    ttesttype = 'indepsamplesT';
    
    cfg_ga = [];
    cfg_ga.keepindividual = 'yes';
    cond1data = ft_timelockgrandaverage(cfg_ga, tldata{1:numsubj1,1});
    cond2data = ft_timelockgrandaverage(cfg_ga, tldata{numsubj1+1:end,1});
    
    if size(tldata,2) > 1
        cond1sub = ft_timelockgrandaverage(cfg_ga, tldata{1:numsubj1,2});
        cond2sub = ft_timelockgrandaverage(cfg_ga, tldata{numsubj1+1:end,2});
        cond1data.individual = cond1data.individual - cond1sub.individual;
        cond2data.individual = cond2data.individual - cond2sub.individual;
    end
    
    cond1data.avg = squeeze(mean(cond1data.individual,1));
    cond2data.avg = squeeze(mean(cond2data.individual,1));
    
    design = zeros(1,numsubj);
    design(1,1:numsubj1) = 1;
    design(1,numsubj1+1:end)= 2;
    cfg.ivar  = 1;                   % number or list with indices, independent variable(s)
end

if (strcmp(statmode,'cond') || strcmp(statmode,'subj')) && (strcmp(param.testcnv,'on') || strcmp(param.testlat,'on'))
    timeidx = cond1data.time >= latency(1)+timeshift & cond1data.time <= latency(2)+timeshift;
    for ind = 1:size(cond1data.individual,1)
        for chan = 1:length(cond1data.label)
            if strcmp(param.testcnv,'on')
                summinfo = polyfit(cond1data.time(timeidx),squeeze(cond1data.individual(ind,chan,timeidx))',1);
            elseif strcmp(param.testlat,'on')
                summinfo = calclat(cond1data.time(timeidx),squeeze(cond1data.individual(ind,chan,timeidx))',50);
            end
            cond1data.individual(ind,chan,1) = summinfo(1);
        end
    end
    cond1data.time = 0;
    cond1data.individual = cond1data.individual(:,:,1);
    cond1data.avg = squeeze(mean(cond1data.individual,1))';
    
    timeidx = cond2data.time >= latency(1)+timeshift & cond2data.time <= latency(2)+timeshift;
    for ind = 1:size(cond2data.individual,1)
        for chan = 1:length(cond2data.label)
            if strcmp(param.testcnv,'on')
                summinfo = polyfit(cond2data.time(timeidx),squeeze(cond2data.individual(ind,chan,timeidx))',1);
            elseif strcmp(param.testlat,'on')
                summinfo = calclat(cond2data.time(timeidx),squeeze(cond2data.individual(ind,chan,timeidx))',50);
            end
            cond2data.individual(ind,chan,1) = summinfo(1);
        end
    end
    cond2data.time = 0;
    cond2data.individual = cond2data.individual(:,:,1);
    cond2data.avg = squeeze(mean(cond2data.individual,1))';
end

cfg.design = design;             % design matrix
cfg.statistic = ttesttype;

diffcond = cond1data;
diffcond.cond1avg = cond1data.avg;
diffcond.cond2avg = cond2data.avg;
diffcond.avg = cond1data.avg - cond2data.avg;


fprintf('\nComparing conditions using %d-tailed %s test\nat alpha of %.2f between %.2f-%.2f sec.\n\n', param.ttesttail, ttesttype, param.alpha, latency);
cfg.latency = latency + timeshift;            % time interval over which the experimental conditions must be compared (in seconds)

cfg.feedback = 'textbar';

[stat] = ft_timelockstatistics(cfg, cond1data, cond2data);
stat.chanlocs = chanlocs;

save2file = sprintf('%s_%s_%s-%s.mat',statmode,num2str(subjinfo),condlist{1},condlist{2});
stat.cfg = cfg;
stat.condlist = condlist;
stat.diffcond = diffcond;
stat.timeshift = timeshift;
stat.statmode = statmode;
stat.subjinfo = subjinfo;

plotclusters(stat);

if nargout == 0
    save(save2file, 'stat');
end

function EEG = convertogfp(EEG)

EEG.nbchan = 1;
EEG.trials = 1;
EEG.chanlocs = EEG.chanlocs(find(strcmp('Cz',{EEG.chanlocs.labels})));
[~, EEG.data] = evalc('eeg_gfp(mean(EEG.data,3)'')''');

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
