function plot_topo(P,D)
% requires:
    % D.DAT (output from plotprepare_)
    % P.WFrows
    % P.Nplots
    % P.Fi_ind
    % P.fi
    % P.cond
    % P.chanlocs
    % P.Eval
    % P.EEGtimes
    
% average ERPs over non-unique rows
ind = unique(D.WFrows,'stable')';
DATa = cell(length(ind),1);
for i = ind
    DATa{i}=mean(cat(3,D.DAT{D.WFrows==i}),3);
end
f1=figure
pln=0;

if isfield(P,'selectlev')
    D.Fi_ind = find(D.Fi(:,1)==P.selectlev); % indices of Fi (and wf) for each plot
%else
%    D.Fi_ind = find(D.Fi(:,1)); % indices of Fi (and wf) for each plot
end

y=DATa(D.Fi_ind);
y=y(D.fi);
if iscell(P.cval)
    cond = P.cval{1}(P.cval{2});
    condind=nan(1,length(D.cond));
    for cn = 1:length(cond) 
        ind = ismember(D.cond,cond{cn});
        condind(ind) = cn;
    end
else
    [cond,~,condind]=unique(D.cond,'stable');
    unicond = unique(condind)';
end

% gather data
plotchans=1:length(D.chanlocs);
plotchans(P.no_plot_ele)=[];
dat=struct;
for cn = 1:length(cond) 
    peakdata=y(condind==cn);
    peakdata = mean(cat(3,peakdata{:}),3);
    %[~,markchans] = intersect(plotchans,tp);
    if any(P.topo_subtimewin) && length(P.topo_subtimewin)==1 % multiple plots within a range
        cluswin = D.E_val(end)-D.E_val(1);
        nlat = floor(cluswin/P.topo_subtimewin);
        if nlat>1
            lats = linspace(D.E_val(1),D.E_val(end),nlat);
            lats = dsearchn(D.EEGtimes',lats');
        else
            lats = [find(D.EEGtimes==D.E_val(1)) find(D.EEGtimes==D.E_val(end))];
        end

        for ln = 1:length(lats)-1
            pln = pln+1;
            lat=lats(ln:ln+1);
            dat(pln).peakdata = mean(peakdata(:,min(lat):max(lat)),2);
            ax(pln) = subplot(length(cond),length(lats)-1,pln); 
        end
    elseif any(P.topo_subtimewin) && length(P.topo_subtimewin)==2 % specified time window
        pln = pln+1;
        lat = dsearchn(D.EEGtimes',P.topo_subtimewin');
        dat(pln).peakdata = mean(peakdata(:,min(lat):max(lat)),2);
        ax(pln) = subplot(length(unicond),1,pln); 
    else
        pln = pln+1;
        lat = find(D.EEGtimes==D.E_val(1));
        dat(pln).peakdata = mean(peakdata(:,min(lat):max(lat)),2);
        ax(pln) = subplot(length(unicond),1,pln); 
    end
    %title(num2str(D.EEGtimes(lat)),'fontsize',P.fontsize)
end

% equate scales
alldat=[dat(:).peakdata];
plotmin = min(alldat(:));
plotmax = max(alldat(:));

plotchans=1:length(D.chanlocs);
plotchans(P.no_plot_ele)=[];
for pln = 1:length(dat) 
    axes(ax(pln))
    plottopotype(dat(pln).peakdata,D.chanlocs,plotchans,P.topotype,[plotmin plotmax]) 
end
linkaxes(ax,'xy')
end

function plottopotype(dat,chanlocs,plotchans,topotype,minmax)

switch topotype
    case 'eeglab'
        topoplot(dat, chanlocs,'maplimits',minmax,'electrodes','off','style','map','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
    case 'pcolor' % UNFINISHED
        topo = nanmean(R(f2).Vdv(:,:,tw),3);
        figure
        set(gcf, 'Position', [400, 400, 350, 300])
        %imagesc(rot90(topo))
        pcolor(rot90(topo,3)), shading interp
        axis off
        colormap(jet)
        title([field{f} ': time of max ' field{f2} ' for image ' num2str(xval_peak) 'ms'])
        colorbar
end
   
%tightfig(f1)
end