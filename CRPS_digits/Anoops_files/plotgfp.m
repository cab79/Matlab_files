function plotgfp(stat,varargin)

loadpaths

if strcmp(stat.condlist{end},'base')
    condlist = stat.condlist(1:end-1);
else
    condlist = stat.condlist(1:2);
end

if ~isfield(stat,'inddata')
    stat.inddata{1} = stat.condavg(:,:,1);
    stat.inddata{2} = stat.condavg(:,:,2);
end

colorlist = {
    'explicit'     [0         0    1.0000]
    'implicit'     [0    0.5000         0]
    'distractor'   [1.0000    0         0]
    'target'      [0    0.7500    0.7500]
    };

param = finputcheck(varargin, { 'ylim', 'real', [], [-5 15]; ...
    'fontsize','integer', [], 26; ...
    'legendstrings', 'cell', {}, condlist; ...
    'legendposition', 'string', {}, 'NorthEast'; ...
    'ttesttail', 'integer', [-1 0 1], 0; ...
    'plotinfo', 'string', {'on','off'}, 'on'; ...
    });

%% figure plotting

fontname = 'Helvetica';
linewidth = 1;

if length(condlist) == 2
    figname = sprintf('%s-%s',stat.condlist{1},stat.condlist{2});
else
    figname = sprintf('%s',stat.condlist{1});
end

figfile = sprintf('figures/%s_%s_%s_gfp',stat.statmode,num2str(stat.subjinfo),figname);

figure('Name',figname,'Color','white','FileName',[figfile '.fig']);
figpos = get(gcf,'Position');
set(gcf,'Position',[50 50 figpos(3) figpos(3)]);

if ~isempty(stat.pclust) && length(stat.pclust.win(1):4:stat.pclust.win(2)) < 10
    stat.pclust = struct([]);
end

if ~isempty(stat.pclust)
    latpnt = find(stat.times-stat.timeshift >= stat.pclust(1).win(1) & stat.times-stat.timeshift <= stat.pclust(end).win(2));
else
    latpnt = find(stat.times-stat.timeshift >= stat.param.latency(1) & stat.times-stat.timeshift <= stat.param.latency(2));
end

%pick time point at max of condition 1
[~, maxidx] = max(stat.condgfp(1,latpnt,1),[],2);

%pick time point at max of difference
%[~, maxidx] = max(stat.gfpdiff(1,latpnt),[],2);

%pick time point at max of t-statistic
%[~, maxidx] = max(stat.valu(latpnt));

%pick time point at min of p-value
%[~, maxidx] = min(stat.pprob(latpnt));

plotpnt = latpnt(1)-1+maxidx;

subplot(2,2,1:2);
plotvals = mean(stat.inddata{1}(:,plotpnt,:),3);
topoplot(plotvals,stat.chanlocs,'electrodes','labels');
cb_h = colorbar('FontSize',param.fontsize);
cb_labels = num2cell(get(cb_h,'YTickLabel'),2);
cb_labels{1} = [cb_labels{1} ' uV'];
set(cb_h,'YTickLabel',cb_labels);

if ~isempty(stat.pclust)
        text(0,-0.7,sprintf('%dms\nclust. t = %.1f, p = %.3f', ...
            round(stat.times(plotpnt)), ...
            stat.pclust(1).tstat, stat.pclust(1).prob),...
            'FontSize',param.fontsize,'FontName',fontname,'HorizontalAlignment','center');
else
        text(0,-0.7,sprintf('%dms', round(stat.times(plotpnt))), ...
            'FontSize',param.fontsize,'FontName',fontname,'HorizontalAlignment','center');
end

subplot(2,2,3:4);
curcolororder = get(gca,'ColorOrder');
colororder = zeros(length(param.legendstrings),3);
for str = 1:length(param.legendstrings)
    cidx = strcmp(param.legendstrings{str},colorlist(:,1));
    if sum(cidx) == 1
        colororder(str,:) = colorlist{cidx,2};
    else
        colororder(str,:) = curcolororder(str,:);
    end
end

% colororder = cat(1,colororder,[0 0 0]);
set(gca,'ColorOrder',colororder);
hold all;

plotdata = squeeze(stat.condgfp(1,:,1:length(condlist)));
if length(condlist) == 2
    plotdata = plotdata';
end

%plotdata = stat.gfpdiff(1,:);
%param.legendstrings{end+1} = 'difference';

% if strcmp(stat.statmode,'trial')
    plot((stat.times(1):1000/stat.srate:stat.times(end))-stat.timeshift,plotdata,'LineWidth',linewidth*1.5);
% elseif strcmp(stat.statmode,'cond') || strcmp(stat.statmode,'subj')
%     H = shadedErrorBar((stat.times(1):1000/stat.srate:stat.times(end))-stat.timeshift,plotdata(1,:),...
%         std(stat.indgfp{1}'),{'LineWidth',linewidth*1.5,'Color',colororder(1,:)});
%     hAnnotation = get(H.patch,'Annotation');
%     hLegendEntry = get(hAnnotation','LegendInformation');
%     set(hLegendEntry,'IconDisplayStyle','off');
%     for e = 1:length(H.edge)
%         hAnnotation = get(H.edge(e),'Annotation');
%         hLegendEntry = get(hAnnotation','LegendInformation');
%         set(hLegendEntry,'IconDisplayStyle','off');
%     end
% 
%     if length(condlist) == 2
%         H = shadedErrorBar((stat.times(1):1000/stat.srate:stat.times(end))-stat.timeshift,plotdata(2,:),...
%             std(stat.indgfp{2}'),{'LineWidth',linewidth*1.5,'Color',colororder(2,:)});
%         hAnnotation = get(H.patch,'Annotation');
%         hLegendEntry = get(hAnnotation','LegendInformation');
%         set(hLegendEntry,'IconDisplayStyle','off');
%         for e = 1:length(H.edge)
%             hAnnotation = get(H.edge(e),'Annotation');
%             hLegendEntry = get(hAnnotation','LegendInformation');
%             set(hLegendEntry,'IconDisplayStyle','off');
%         end
%     end
% end

set(gca,'XLim',[stat.times(1) stat.times(end)]-stat.timeshift,'XTick',stat.times(1)-stat.timeshift:200:stat.times(end)-stat.timeshift,'YLim',param.ylim,...
    'FontSize',param.fontsize,'FontName',fontname);
line([0 0],param.ylim,'Color','black','LineStyle',':','LineWidth',linewidth);
line([stat.times(1) stat.times(end)]-stat.timeshift,[0 0],'Color','black','LineStyle',':','LineWidth',linewidth);
line([stat.times(plotpnt) stat.times(plotpnt)]-stat.timeshift,param.ylim,'Color','red','LineWidth',linewidth,'LineStyle','--');
if strcmp(param.plotinfo,'on')
    xlabel('Time (ms) ','FontSize',param.fontsize,'FontName',fontname);
    ylabel('Global field power ','FontSize',param.fontsize,'FontName',fontname);
else
    xlabel(' ','FontSize',param.fontsize,'FontName',fontname);
    ylabel(' ','FontSize',param.fontsize,'FontName',fontname);
end
legend(param.legendstrings,'Location',param.legendposition);

%% plot clusters

if param.ttesttail >= 0
    if ~isempty(stat.pclust)
        for p = 1:length(stat.pclust)
            line([stat.pclust(p).win(1) stat.pclust(p).win(2)],[0 0],'Color','blue','LineWidth',8);
        end
    else
        fprintf('No positive clusters found.\n');
    end
end

if param.ttesttail <= 0
    if ~isempty(stat.nclust)
        for n = 1:length(stat.nclust)
            line([stat.nclust(n).win(1) stat.nclust(n).win(2)],[0 0],'Color','red','LineWidth',8);
        end
    else
        fprintf('No negative clusters found.\n');
    end
    
end
    
set(gcf,'Color','white');
%saveas(gcf,[figfile '.fig']);
% export_fig(gcf,[figfile '.eps']);
