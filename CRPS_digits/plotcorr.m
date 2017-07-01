function plotcorr(stat,data,cov)

posclustidx = [];
if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
    for cidx = 1:length(stat.posclusters)
        if stat.posclusters(cidx).prob < stat.cfg.alpha && isempty(posclustidx) ...
                || (~isempty(posclustidx) && stat.posclusters(cidx).prob < stat.posclusters(posclustidx).prob)
            posclustidx = cidx;
        end
    end
end

negclustidx = [];
if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
    for cidx = 1:length(stat.negclusters)
        if stat.negclusters(cidx).prob < stat.cfg.alpha && isempty(negclustidx) ...
                || (~isempty(negclustidx) && stat.negclusters(cidx).prob < stat.negclusters(negclustidx).prob)
            negclustidx = cidx;
        end
    end
end

figure('Name','scatter','Color','white');
if ~isempty(posclustidx)
    win = find(ismember(data.time,stat.time(find(stat.posclusterslabelmat))));
    y=squeeze(mean(data.individual(:,:,win),3))
    x=cov;
    scatter(x,y);
    lsline;
end