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
DATa = cell(length(unique(D.WFrows)),1);
for i = unique(D.WFrows)'
    DATa{i}=mean(cat(3,D.DAT{D.WFrows==i}),3);
end
f1=figure
pln=0;

y=DATa(D.Fi_ind);
y=y(D.fi);
[~,~,condind]=unique(D.cond,'stable');
unicond = unique(condind)';
for cn = unicond
    peakdata=y(condind==cn);
    peakdata = mean(cat(3,peakdata{:}),3);
    plotchans=1:length(D.chanlocs);
    plotchans(P.no_plot_ele)=[];
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
            subplot(length(unicond),length(lats)-1,pln); topoplot(mean(peakdata(:,lat(1):lat(2)),2), D.chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
            title(num2str(EEGtimes(lat)'))
        end
    elseif any(P.topo_subtimewin) && length(P.topo_subtimewin)==2 % specified time window
        pln = pln+1;
        lat = dsearchn(D.EEGtimes',P.topo_subtimewin');
        subplot(length(unicond),1,pln); topoplot(mean(peakdata(:,lat),2), D.chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
        title(num2str(D.EEGtimes(lat)'))
    else
        pln = pln+1;
        lat = find(D.EEGtimes==D.E_val(1));
        subplot(length(unicond),1,pln); topoplot(mean(peakdata(:,lat),2), D.chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
        title(num2str(D.EEGtimes(lat)'))
    end
end
tightfig(f1)