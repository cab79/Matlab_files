function CORE_eeg_plot_stats(name,xticks,topo_range,t,fdr_t,chanlocs)
figure('name',name)
trange = dsearchn(xticks',[topo_range(1);topo_range(2)]);
clim=[min(t(:)) max(t(:))];
hold on
subplot(2,2,1)
imagesc(xticks,[],t,clim); 
[~,mi]=max(mean(abs(t(:,trange(1):trange(2))),1));
mi=mi+trange(1)-1;
line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
title(['t values for two-sample t-test']);
colorbar
subplot(2,2,2)
topoplot(t(:,mi),chanlocs,'maplimits',clim);
subplot(2,2,3)
clim=[min(fdr_t(:)) max(fdr_t(:))];
if diff(clim)==0; clim(2) = clim(2)+abs(clim(2)*1.01); end
if ~any(clim)
    clim=[-5 5];
end
imagesc(xticks,[],fdr_t,clim); 
[~,mi]=max(mean(abs(fdr_t(:,trange(1):trange(2))),1));
mi=mi+trange(1)-1;
line(xticks(mi*[1 1]),[1 92],'color','k','linewidth',2)
title(['fdr thresholded t values']);
colorbar
subplot(2,2,4)
topoplot(fdr_t(:,mi),chanlocs,'maplimits',clim);
