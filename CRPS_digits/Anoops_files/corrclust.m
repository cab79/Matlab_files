function stat = corrclust(stat,varargin)

param = finputcheck(varargin, {
    'alpha' , 'real' , [], 0.05; ...
    'ttesttail' , 'integer' , [-1 0 1], 0; ...
    'clustsize', 'integer', [], 1; ...
    });

if param.ttesttail == 0
    param.alpha = param.alpha / 2;
end

valu = zeros(1,size(stat.gfpdiff,2));
pprob = ones(1,size(stat.gfpdiff,2));
nprob = ones(1,size(stat.gfpdiff,2));
pmask = zeros(size(pprob));
nmask = zeros(size(nprob));
stat.pclust = struct([]);
stat.nclust = struct([]);

corrwin = find(stat.times >= stat.param.latency(1) & stat.times <= stat.param.latency(2));
stat.pdist = max(stat.gfpdiff(:,corrwin),[],2);
stat.ndist = min(stat.gfpdiff(:,corrwin),[],2);

%% identfy clusters
h_wait = waitbar(0,'Please wait...');
set(h_wait,'Name',[mfilename ' progress']);

for n = 1:size(stat.gfpdiff,1)
    if n > 1
        waitbar(n/size(stat.gfpdiff,1),h_wait,sprintf('Permutation %d...',n-1));
    end
    valu(:) = 0;
    pprob(:) = 1;
    nprob(:) = 1;
    pmask(:) = 0;
    nmask(:) = 0;
    for p = corrwin
        valu(p) = (stat.gfpdiff(n,p) - mean(stat.gfpdiff(:,p)))/...
            (std(stat.gfpdiff(:,p))/sqrt(size(stat.gfpdiff,1)));
        pprob(p) = sum(stat.pdist >= stat.gfpdiff(n,p))/size(stat.gfpdiff,1);
        nprob(p) = sum(stat.ndist <= stat.gfpdiff(n,p))/size(stat.gfpdiff,1);
    end
    
    pmask(pprob < param.alpha) = 1;
    nmask(nprob < param.alpha) = 1;
    
    if n == 1
        stat.valu = valu;
        stat.pprob = pprob;
        stat.nprob = nprob;
        stat.pmask = pmask;
        stat.nmask = nmask;
    end
    
    stat.pclust(n).tstat = 0;
    stat.nclust(n).tstat = 0;
    for p = 2:length(stat.times)
        if pmask(p) == 1 && pmask(p-1) == 0
            pstart = p;
        elseif (pmask(p) == 0 || p == length(stat.times)) && pmask(p-1) == 1
            pend = p;
            tstat = sum(valu(pstart:pend-1));
            if tstat > stat.pclust(n).tstat && length(pstart:pend-1) >= param.clustsize
                stat.pclust(n).tstat = tstat;
                stat.pclust(n).win = [stat.times(pstart) stat.times(pend-1)]-stat.timeshift;
            end
        end
        
        if nmask(p) == 1 && nmask(p-1) == 0
            nstart = p;
        elseif (nmask(p) == 0 || p == length(stat.times)) && nmask(p-1) == 1
            nend = p;
            tstat = sum(valu(nstart:nend-1));
            if tstat > stat.nclust(n).tstat && length(nstart:nend-1) >= param.clustsize
                stat.nclust(n).tstat = tstat;
                stat.nclust(n).win = [stat.times(nstart) stat.times(nend-1)]-stat.timeshift;
            end
        end
    end
end
close(h_wait);

randtstat = cell2mat({stat.pclust(:).tstat});
stat.pclust(1).prob = sum(randtstat >= stat.pclust(1).tstat)/length(randtstat);
stat.pclust(1).tstat = (stat.pclust(1).tstat - mean(randtstat)) / (std(randtstat)/sqrt(length(randtstat)));
stat.pclust = stat.pclust(1);
if stat.pclust.prob > param.alpha
    stat.pclust = struct([]);
end

randtstat = cell2mat({stat.nclust(:).tstat});
stat.nclust(1).prob = sum(randtstat <= stat.nclust(1).tstat)/length(randtstat);
stat.nclust(1).tstat = (stat.nclust(1).tstat - mean(randtstat)) / (std(randtstat)/sqrt(length(randtstat)));
stat.nclust = stat.nclust(1);
if stat.nclust.prob > param.alpha
    stat.nclust = struct([]);
end

paramlist = fieldnames(param);
for p = 1:length(paramlist)
    stat.param.(paramlist{p}) = param.(paramlist{p});
end
