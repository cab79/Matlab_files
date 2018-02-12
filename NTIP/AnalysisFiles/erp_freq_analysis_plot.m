function S=erp_freq_analysis_plot(S)

% GET FILE LIST
S.filepath = S.freqpath;
S.loadext='Freq.mat';
S = getfilelist(S);

% select session, block, frequency
S.fn = dsearchn(S.freqsrange',S.freqselect);

if ~exist('chanlocs','var')
    load(fullfile(S.setpath,'chanlocs.mat'));
    S.chanlocs = chanlocs;
end
S.inclchan = 1:length(S.chanlocs);
S.inclchan(S.exclchan) = [];
for s = 1:length(S.subjects)
    %if participant_eventorder(s)==-1
    %    event_order = fliplr(event_reorder);
    %else
    %    event_order = event_reorder;
    %end
    
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,S.subjects{s})))));
    % FIND THE FILES FOR THIS BLOCK
    blockfiles = subfiles(find(not(cellfun('isempty', strfind(subfiles,S.blocks{S.blockselect})))));
    
    % get the data
    for c = 1:length(S.conds)
        % get filename
        file = blockfiles(find(not(cellfun('isempty', strfind(blockfiles,S.conds{c})))));
        if length(file)~=1
            error('file name not uniquely specified')
        end
        % load data
        load(fullfile(S.filepath,file{:}));
        % average data over selected conditions
        if isstruct(fdata)
            dat = fdata.powspctrm; %trials,chans,freqs
            madat = fdata.madata;
            act_freq = fdata.freq;
        elseif iscell(fdata)
            clear datc madatc
            for e = 1:length(S.eventselect)
                datc{e} = fdata{S.eventselect(e)}.powspctrm;
                madatc{e} = fdata{S.eventselect(e)}.madata;
            end
            dat = cat(1,datc{:});
            madat = cat(4,madatc{:});
            act_freq = fdata{1}.freq;
        end
        S.meandat{c,1} = squeeze(mean(dat,1));
        S.movmeandat{c,1} = squeeze(mean(madat,4));
    end

    if ~all(ismember(S.freqsrange,act_freq))
        error('output frequencies do not match the settings specified')
    end


    S.plotdat=S.meandat;
    S.plotmovdat=S.movmeandat;
    
    [condlabels,~]=strtok(S.conds,'_');

    % plot topographies for each condition
    ttl = condlabels;
    lgnd = {};
    scalezero = 0;
    plot_freq(S,ttl,lgnd,scalezero)
    % all frequencies for each condition
    ttl = 'Power at each frequency (no baseline correction)';
    lgnd = condlabels;
    plot_allfreq(S,ttl,lgnd)
    % over time for each condition
    ttl = 'Power over time (no baseline correction)';
    lgnd = condlabels;
    plot_movavg(S,ttl,lgnd,0)

    % Ratio baseline correction
    S=event_operation(S,1);
    % Ratio baseline correction: plot topographies for each condition
    ttl = condlabels(1:5);
    lgnd = {};
    scalezero = 0;
    plot_freq(S,ttl,lgnd,scalezero)
    % Ratio baseline correction: plot all frequencies for each condition
    ttl = 'Baseline corrected: Power at each frequency';
    lgnd = condlabels(1:4);
    plot_allfreq(S,ttl,lgnd)
    % Ratio baseline correction: over time for each condition
    ttl = 'Baseline corrected: Power over time';
    lgnd = condlabels(1:4);
    plot_movavg(S,ttl,lgnd,0)

    % condition subtraction
    S=event_operation(S,2);
    % plot topographies for each condition
    ttl = {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
    lgnd = {};
    scalezero = 0;
    plot_freq(S,ttl,lgnd,scalezero)
    % Ratio baseline correction: plot all frequencies for each condition
    ttl = 'Subtracted: Entrainment effect at each frequency';
    lgnd = {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
    plot_allfreq(S,ttl,lgnd)
    % Ratio baseline correction: over time for each condition
    ttl = 'Subtracted: Entrainment effect over time';
    lgnd = {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
    plot_movavg(S,ttl,lgnd,0)

    %S=freq_operation(S,3);
    %plot_freq(S,0)
    %plot_movavg(S,{'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball'})


end

function S=freq_operation(S,op)
opdat = S.plotdat; 
opmovdat = S.plotmovdat;
S.plotdat = [];
S.plotmovdat = [];
for c = 1:length(opdat)
    for fi = 1:length(S.op(op).grouping)
        fni = 1:length(S.op(op).grouping)
        fni(fi)=[];
        S.plotdat{c}(:,fi) = operation(opdat{c}(:,fi),mean(opdat{c}(:,fni),2),S.op(op).operation);
        if S.mov_avg_trials
            S.plotmovdat{c}(:,fi,:) = operation(opmovdat{c}(:,fi,:),mean(opmovdat{c}(:,fni,:),2),S.op(op).operation);
        end
    end
end

function S=event_operation(S,op)
opdat = S.plotdat; 
opmovdat = S.plotmovdat;
S.plotdat = [];
S.plotmovdat = []; 
for su = 1:length(S.op(op).grouping)
    nFreq = size(opdat{S.op(op).grouping(su,1)},2);
    dat1 = opdat{S.op(op).grouping(su,1)};
    dat2 = opdat{S.op(op).grouping(su,2)};
    S.plotdat{su} = operation(dat1,dat2,S.op(op).operation);
    if S.mov_avg_trials
        dat1 = opmovdat{S.op(op).grouping(su,1)};
        dat2 = repmat(mean(opmovdat{S.op(op).grouping(su,2)},3),1,1,size(dat1,3));
        S.plotmovdat{su}= operation(dat1,dat2,S.op(op).operation);
    end
end

function data=operation(data,meanVals,baselinetype)
if (strcmp(baselinetype, 'subtract'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
  data = (data - meanVals) ./ (data + meanVals);
elseif (strcmp(baselinetype, 'db'))
  data = 10*log10(data ./ meanVals);
else
  error('unsupported method for baseline normalization: %s', baselinetype);
end

function plot_freq(S,ttl,lgnd,scalezero)
% plot a frequency
figure
datmat = cat(3,S.plotdat{:});
maxdat = max(abs(datmat(:)));
minmax = [-maxdat,maxdat]+scalezero;
gsf1 = 0;
for c = 1:length(S.plotdat)
    gsf1=gsf1+1;
    dat = S.plotdat{c}(:,S.fn);
    subplot(1,length(S.plotdat),gsf1)
    topoplot(dat(S.inclchan),S.chanlocs(S.inclchan),'maplimits',minmax);
    cbar('vert',0,minmax);
    title(ttl{c});
end

function plot_allfreq(S,ttl,lgnd)
% plot all frequencies at maximal electrodes
% find electrodes with maximal power at S.analysistypeed frequency
nele = 5;
datmat = cat(3,S.plotdat{:});
dat = mean(datmat(:,S.fn),3);
[~,maxele] = sort(dat,'descend');
maxele=maxele(1:nele);
%plot
fplot=[];
for c = 1:length(S.plotdat)
    for fi = 1:length(S.freqsrange)
        fplot(fi,c) = mean(S.plotdat{c}(maxele,fi),1);
    end
end
figure
plot(S.freqsrange,fplot);
legend(lgnd)
title(ttl)
xlabel('Frequency')

function plot_movavg(S,ttl,lgnd,eq)
% plot trends for selected freq only
nele = 5;
datmat = cat(3,S.plotdat{:});
dat = mean(datmat(:,S.fn),3);
[~,maxele] = sort(dat,'descend');
maxele=maxele(1:nele);
cplot={};
minlen = inf;
maxlen = 0;
for c = 1:length(S.plotmovdat)
    cplot{c} = squeeze(mean(S.plotmovdat{c}(maxele,S.fn,:),1));
    minlen = min(minlen,length(cplot{c}));
    maxlen = max(maxlen,length(cplot{c}));
end
if eq
    % equate lengths
    for c = 1:length(S.plotmovdat)
        cplot{c} = cplot{c}(1:minlen);
    end
else % pad
    for c = 1:length(S.plotmovdat)
        if length(cplot{c})<maxlen
            cplot{c}(end+1:maxlen) = nan;
        end
    end
end
fplot = cat(2,cplot{:});
figure
plot(1:size(fplot,1),fplot);
legend(lgnd)
title(ttl)
xlabel(['Time (N x ' num2str(S.mov_avg_trials_step) ' trials)'])
