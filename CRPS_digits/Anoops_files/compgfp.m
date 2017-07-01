function stat = compgfp(subjinfo,condlist,varargin)

loadpaths

global chanidx

timeshift = 0; %milliseconds

param = finputcheck(varargin, {
    'numrand', 'integer', [], 1000; ...
    'latency', 'real', [], []; ...
    'chanlist', 'cell', {}, {}; ...
    'wori', 'cell', {}, cell(1,length(condlist)), ...
    });

if isempty(param.latency)
    param.latency = [0 EEG.times(end)-timeshift];
end

%% SELECTION OF SUBJECTS AND LOADING OF DATA
loadsubj

if ischar(condlist)
    condlist = {condlist,'base'};
elseif iscell(condlist) && length(condlist) == 1
    condlist{2} = 'base';
end

if ischar(subjinfo)
    %%%% perform single-trial statistics
    subjlist = {subjinfo};
    subjcond = condlist;
    statmode = 'trial';
    
elseif isnumeric(subjinfo) && length(subjinfo) == 1
    %%%% perform within-subject statistics
    subjlist1 = subjlists{subjinfo(1)};
    subjlist2 = subjlists{subjinfo(2)};
    subjlist = cat(2,subjlist1,subjlist2);
    if length(condlist) == 3
        condlist = {sprintf('%s-%s',condlist{1},condlist{3}),sprintf('%s-%s',condlist{2},condlist{3})};
    end
    statmode = 'cond';
    
elseif isnumeric(subjinfo) && (length(subjinfo) == 2 || length(subjinfo) == 4)
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
    statmode = 'subj';
end

numsubj = length(subjlist);
numcond = size(subjcond,2);

conddata = cell(numsubj,numcond);

%% load and prepare individual subject datasets

for s = 1:numsubj
    
    for c = 1:numcond
        EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,c}), 'filepath', filepath);
        EEG = sortchan(EEG);
        
        if s == 1 && c == 1
            chanlocs = EEG.chanlocs;
            times = EEG.times - timeshift;
            corrwin = find(times >= param.latency(1) & times <= param.latency(2));
        end
        % rereference
        %     EEG = rereference(EEG,1);
        
        %%%%% baseline correction relative to whole epoch
        %     EEG = pop_rmbase(EEG,[]);
        %%%%%
        
        conddata{s,c} = EEG;
        
        %         if strcmp(statmode,'trial') && c == numcond
        %             if conddata{s,1}.trials > conddata{s,2}.trials
        %                 fprintf('Equalising trials in condition %s.\n',subjcond{s,1});
        %                 randtrials = 1:conddata{s,1}.trials;%randperm(conddata{s,1}.trials);
        %                 conddata{s,1} = pop_select(conddata{s,1},'trial',randtrials(1:conddata{s,2}.trials));
        %             elseif conddata{s,2}.trials > conddata{s,1}.trials
        %                 fprintf('Equalising trials in condition %s.\n',subjcond{s,2});
        %                 randtrials = 1:conddata{s,2}.trials;%randperm(conddata{s,2}.trials);
        %                 conddata{s,2} = pop_select(conddata{s,2},'trial',randtrials(1:conddata{s,1}.trials));
        %             end
        %         end
        
    end
end

% chanidx = {'E54','E55','E79','E78','E62','E61','E31','E80','E87','E86','E85','E77','E72','E67','E60','E53','E37'};
% 
% for c = 1:length(chanidx)
%     chanidx{c} = find(strcmp(chanidx{c},{chanlocs.labels}));
% end
% chanidx = cell2mat(chanidx);

if strcmp(statmode,'trial')
    inddata{1} = conddata{1}.data;
    inddata{2} = conddata{2}.data;
    indgfp{1} = calcgfp(mean(inddata{1},3),EEG.times);
    indgfp{2} = calcgfp(mean(inddata{2},3),EEG.times);
    mergedata = cat(3,conddata{1,1}.data,conddata{1,2}.data);
    diffcond = mean(inddata{1},3) - mean(inddata{2},3);
    
elseif strcmp(statmode,'cond')
    inddata{1} = zeros(conddata{1,1}.nbchan,conddata{1,1}.pnts,numsubj);
    inddata{2} = zeros(conddata{1,2}.nbchan,conddata{1,2}.pnts,numsubj);
    indgfp{1} = zeros(conddata{1,1}.pnts,numsubj);
    indgfp{2} = zeros(conddata{1,2}.pnts,numsubj);
    
    for s = 1:numsubj
        inddata{1}(:,:,s) = mean(conddata{s,1}.data,3);
        inddata{2}(:,:,s) = mean(conddata{s,2}.data,3);
        
        if size(conddata,2) > 2
            condsub = mean(conddata{s,3}.data,3);
            inddata{1}(:,:,s) = inddata{1}(:,:,s) - condsub;
            inddata{2}(:,:,s) = inddata{2}(:,:,s) - condsub;
        end
        indgfp{1}(:,s) = calcgfp(inddata{1}(:,:,s),EEG.times);
        indgfp{2}(:,s) = calcgfp(inddata{2}(:,:,s),EEG.times);
    end
    mergedata = cat(3,inddata{1},inddata{2});
    diffcond = mean(inddata{1},3) - mean(inddata{2},3);
    
elseif strcmp(statmode,'subj')
    inddata{1} = zeros(conddata{1,1}.nbchan,conddata{1,1}.pnts,numsubj1);
    indgfp{1} = zeros(conddata{1,1}.pnts,numsubj1);
    for s = 1:numsubj1
        inddata{1}(:,:,s) = mean(conddata{s,1}.data,3);
        if size(conddata,2) > 1
            cond1sub = mean(conddata{s,2}.data,3);
            inddata{1}(:,:,s) = inddata{1}(:,:,s) - cond1sub;
        end
        indgfp{1}(:,s) = calcgfp(inddata{1}(:,:,s),EEG.times);
    end
    
    inddata{2} = zeros(conddata{numsubj1+1,1}.nbchan,conddata{numsubj1+1,1}.pnts,numsubj2);
    indgfp{2} = zeros(conddata{numsubj1+1,1}.pnts,numsubj2);
    for s = 1:numsubj2
        inddata{2}(:,:,s) = mean(conddata{numsubj1+s,1}.data,3);
        if size(conddata,2) > 1
            cond2sub = mean(conddata{numsubj1+s,2}.data,3);
            inddata{2}(:,:,s) = inddata{2}(:,:,s) - cond2sub;
        end
        indgfp{2}(:,s) = calcgfp(inddata{2}(:,:,s),EEG.times);
    end
    mergedata = cat(3,inddata{1},inddata{2});
    diffcond = mean(inddata{1},3) - mean(inddata{2},3);
end

gfpdiff = zeros(param.numrand+1,conddata{1,1}.pnts);
stat.condgfp = zeros(param.numrand+1,conddata{1,1}.pnts,numcond);
stat.inddata = inddata;
stat.indgfp = indgfp;

h_wait = waitbar(0,'Please wait...');
set(h_wait,'Name',[mfilename ' progress']);

for n = 1:param.numrand+1
    if n > 1
        waitbar((n-1)/param.numrand,h_wait,sprintf('Permutation %d...',n-1));
        mergedata = mergedata(:,:,randperm(size(mergedata,3)));
        inddata{1} = mergedata(:,:,1:size(inddata{1},3));
        inddata{2} = mergedata(:,:,size(inddata{1},3)+1:end);
    end
    
    cond1gfp = calcgfp(mean(inddata{1},3),EEG.times);
    cond2gfp = calcgfp(mean(inddata{2},3),EEG.times);
    gfpdiff(n,:) = cond1gfp - cond2gfp;
    stat.condgfp(n,:,1) = cond1gfp;
    stat.condgfp(n,:,2) = cond2gfp;
end
close(h_wait);

stat.valu = zeros(1,size(gfpdiff,2));
stat.pprob = ones(1,size(gfpdiff,2));
stat.nprob = ones(1,size(gfpdiff,2));
stat.pdist = max(gfpdiff(2:end,corrwin),[],2);
stat.ndist = min(gfpdiff(2:end,corrwin),[],2);

for p = corrwin
    stat.valu(p) = (gfpdiff(1,p) - mean(gfpdiff(2:end,p)))/(std(gfpdiff(2:end,p))/sqrt(size(gfpdiff,1)-1));
    stat.pprob(p) = sum(stat.pdist >= gfpdiff(1,p))/param.numrand;
    stat.nprob(p) = sum(stat.ndist <= gfpdiff(1,p))/param.numrand;
end

stat.gfpdiff = gfpdiff;
stat.diffcond = diffcond;
stat.times = EEG.times;
stat.condlist = condlist;
stat.timeshift = timeshift;
stat.subjinfo = subjinfo;
stat.statmode = statmode;
stat.param = param;
stat.chanlocs = chanlocs;
stat.srate = EEG.srate;

stat = corrclust(stat);

if nargout == 0
    save2file = sprintf('%s/%s_%s_%s-%s_%d-%d_gfp.mat',filepath,statmode,num2str(subjinfo),...
        condlist{1},condlist{2},param.latency(1),param.latency(2));
    save(save2file,'stat');
end

clear global chanidx
end

function gfp = calcgfp(data,times)
%data in channels x timepoints

global chanidx
if ~isempty(chanidx)
    data = data(chanidx,:);
end
[~,gfp] = evalc('eeg_gfp(data'',0)''');
gfp = rmbase(gfp,[],1:find(times == 0));
end