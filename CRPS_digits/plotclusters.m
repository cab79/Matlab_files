function plotclusters(stat,varargin)

colorlist = {
    'monaural loc. std.'    [0         0    1.0000]
    'monaural loc. dev.'     [0    0.5000         0]
    'monaural glo. std.'   [1.0000    0         0]
    'monaural glo. dev.'    [0    0.7500    0.7500]
    'interaural dev.'  [0.7500    0    0.7500]
    'interaural ctrl.' [0.7500    0.7500    0]
    'attend tones'      [0    0.5000    0.5000]
    'attend sequences'  [0.5000    0    0.5000]
    'interference'      [0    0.2500    0.7500]
    'early glo. std.'   [0.5000    0.5000    0]
    'late glo. std.'    [0.2500    0.5000    0]
    };

param = finputcheck(varargin, {
    'legendstrings', 'cell', {}, stat.condlist; ...
    'legendposition', 'string', {}, 'SouthEast'; ...
    'ylim', 'real', [], [-1.5 1.5]; ...
    'title','string', {}, ' '; ...
    });


fontname = 'Helvetica';
fontsize = 20;
linewidth = 1;

%% plot significant clusters

posclustidxall = [];
if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
    for cidx = 1:length(stat.posclusters)
        if stat.posclusters(cidx).prob < stat.cfg.alpha %&& isempty(posclustidxall) %...
                %|| (~isempty(posclustidxall) && stat.posclusters(cidx).prob < stat.posclusters(posclustidxall).prob)
            posclustidxall = [posclustidxall cidx];
        end
    end
end

negclustidxall = [];
if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
    for cidx = 1:length(stat.negclusters)
        if stat.negclusters(cidx).prob < stat.cfg.alpha %&& isempty(negclustidxall)% ...
                %|| (~isempty(negclustidxall) && stat.negclusters(cidx).prob < stat.negclusters(negclustidxall).prob)
            negclustidxall = [negclustidxall cidx];
        end
    end
end

if stat.cfg.tail >= 0
    for ci = 1:length(posclustidxall)
        posclustidx = posclustidxall(ci);
        fprintf('Plotting positive clusters.\n');

        figfile = sprintf('figures/%s_%s_%s-%s_pos',stat.statmode,num2str(stat.subjinfo),stat.condlist{1},stat.condlist{2});
        figure('Name',sprintf('%s-%s_Positive_Clusters',stat.condlist{1},stat.condlist{2}),'Color','white','FileName',[figfile '.fig']);
        figpos = get(gcf,'Position');
        figpos(1:2) = 20;
        figpos(4) = figpos(3);
        set(gcf,'Position',figpos);

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

        
        win = find(stat.posclusterslabelmat==posclustidx);
        [maxvalT,maxidxT] = max(stat.stat(win));
        tol = 0.0001;
        maxtime = find(abs(stat.time(maxidxT+win(1)-1)-stat.diffcond.time) < tol);
        if size(stat.diffcond.avg,2)>1
            maxval = stat.diffcond.avg(:,maxtime);
            [~,maxchan] = max(maxval);
        else
            maxval = stat.diffcond.avg(maxtime);
            maxchan = 1;
        end
        maxmaxidx = maxidxT;
        
        %clust_t = stat.diffcond.avg(:,find(min(abs(stat.diffcond.time-stat.cfg.latency(1))) == abs(stat.diffcond.time-stat.cfg.latency(1))):...
        %    find(min(abs(stat.diffcond.time-stat.cfg.latency(2))) == abs(stat.diffcond.time-stat.cfg.latency(2))));
        %if ~isempty(posclustidx)
        %    clust_t(~(stat.posclusterslabelmat == posclustidx)) = 0;
        %end
        %[maxval,maxidx] = max(clust_t);
        %[~,maxmaxidx] = max(maxval);
        %maxchan = maxidx(maxmaxidx);
        %maxtime = find(stat.time(maxmaxidx) == stat.diffcond.time);
        

        plotvals = maxval;
        if length(plotvals)==length(stat.chanlocs)
            sp=2;
            subplot(sp,1,1);

            if ~isempty(posclustidx)

                if size(stat.posclusterslabelmat,1)==1; 
                    mask = repmat(stat.posclusterslabelmat(:,maxmaxidx)==posclustidx,1,length(stat.chanlocs));
                else
                    mask = stat.posclusterslabelmat(:,maxmaxidx)==posclustidx;
                end
                topoplot(plotvals,stat.chanlocs, 'maplimits', 'absmax', 'electrodes','off', 'emarker2',{maxchan,'o','green',14,1},...
                    'pmask',mask);
            else
                topoplot(plotvals,stat.chanlocs, 'maplimits', 'absmax', 'electrodes','off', 'emarker2',{maxchan,'o','green',14,1},...
                    'numcontour',0);
            end

            colorbar('FontName',fontname,'FontSize',fontsize);
            title(param.title,'FontName',fontname,'FontSize',fontsize);
        else
            sp=1;
        end

        subplot(sp,1,sp);
        set(gca,'ColorOrder',cat(1,colororder,[0 0 0]));
        hold all;
        
        if size(stat.diffcond.avg,2)>1
            plot(stat.diffcond.time-stat.timeshift,[stat.diffcond.cond1avg(maxchan,:); stat.diffcond.cond2avg(maxchan,:)]','LineWidth',linewidth*1.5);
        else
            plot(stat.diffcond.time-stat.timeshift,[stat.diffcond.cond1avg, stat.diffcond.cond2avg]','LineWidth',linewidth*1.5);
        end
        %ylim = get(gca,'YLim');
        %ylim = ylim*2;
        ylim = param.ylim;
        set(gca,'YLim',ylim,'XLim',[stat.diffcond.time(1) stat.diffcond.time(end)]-stat.timeshift,'XTick',stat.diffcond.time(1)-stat.timeshift:0.2:stat.diffcond.time(end)-stat.timeshift,...
            'FontName',fontname,'FontSize',fontsize);
        legend(param.legendstrings,'Location',param.legendposition);

        line([stat.diffcond.time(1) stat.diffcond.time(end)]-stat.timeshift,[0 0],'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.60 -0.60],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.45 -0.45],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.30 -0.30],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.15 -0.15],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([    0     0],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([stat.diffcond.time(maxtime) stat.diffcond.time(maxtime)]-stat.timeshift,ylim,'LineWidth',linewidth,'LineStyle','--','Color','red');
        xlabel('Time (s) ','FontName',fontname,'FontSize',fontsize);
        ylabel('Amplitude (uV)','FontName',fontname,'FontSize',fontsize);
        box off
        if ~isempty(posclustidx)
            clustwinidx = find(maxval);
            %         rectangle('Position',[stat.time(clustwinidx(1))-stat.timeshift ylim(1) ...
            %             stat.time(clustwinidx(end))-stat.time(clustwinidx(1)) ylim(2)-ylim(1)],'LineStyle','--','EdgeColor','black','LineWidth',linewidth);
            line([stat.time(clustwinidx(1)) stat.time(clustwinidx(end))]-stat.timeshift,[0 0],'Color','blue','LineWidth',8);
            title(sprintf('Cluster @ %.3f sec (t = %.2f, p = %.3f)', stat.diffcond.time(maxtime)-stat.timeshift,stat.posclusters(posclustidx).clusterstat,stat.posclusters(posclustidx).prob),...
                'FontName',fontname,'FontSize',fontsize);
            line([stat.time(win(1)) stat.time(win(end))]-stat.timeshift,[0 0],'Color','red','LineWidth',8);
        else
            title(sprintf('%.3f sec', stat.diffcond.time(maxtime)-stat.timeshift),'FontName',fontname,'FontSize',fontsize);
        end

        set(gcf,'Color','white');
    %     export_fig(gcf,[figfile '.eps']);
    end
else
    fprintf('No significant positive clusters found.\n');
end

if stat.cfg.tail <= 0
    for ci = 1:length(negclustidxall)
        negclustidx = negclustidxall(ci);
        fprintf('Plotting negative clusters.\n');

        figfile = sprintf('figures/%s_%s_%s-%s_neg',stat.statmode,num2str(stat.subjinfo),stat.condlist{1},stat.condlist{2});
        figure('Name',sprintf('%s-%s_Negative_Clusters',stat.condlist{1},stat.condlist{2}),'Color','white','FileName',[figfile '.fig']);
        figpos = get(gcf,'Position');
        figpos(1:2) = 20;
        figpos(4) = figpos(3);
        set(gcf,'Position',figpos);

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

        set(gca,'ColorOrder',cat(1,colororder,[0 0 0]));
        hold all;

        win = find(stat.negclusterslabelmat==negclustidx);
        [maxvalT,maxidxT] = max(stat.stat(win));
        tol = 0.0001;
        maxtime = find(abs(stat.time(maxidxT+win(1)-1)-stat.diffcond.time) < tol);
        if size(stat.diffcond.avg,2)>1
            maxval = stat.diffcond.avg(:,maxtime);
            [~,maxchan] = max(maxval);
        else
            maxval = stat.diffcond.avg(maxtime);
            maxchan = 1;
        end
        maxmaxidx = maxidxT;
        
        %clust_t = stat.diffcond.avg(:,find(min(abs(stat.diffcond.time-stat.cfg.latency(1))) == abs(stat.diffcond.time-stat.cfg.latency(1))):...
        %    find(min(abs(stat.diffcond.time-stat.cfg.latency(2))) == abs(stat.diffcond.time-stat.cfg.latency(2))));
        %if ~isempty(negclustidx)
        %    clust_t(~(stat.negclusterslabelmat == negclustidx)) = 0;
        %end
        %[maxval,maxidx] = max(clust_t);
        %[~,maxmaxidx] = max(maxval);
        %maxchan = maxidx(maxmaxidx);
        %maxtime = find(stat.time(maxmaxidx) == stat.diffcond.time);





        plotvals = maxval;
        if length(plotvals)==length(stat.chanlocs)
            sp=2;
            subplot(sp,1,1);

            %plot cluster with mask
            if ~isempty(negclustidx)
                if size(stat.negclusterslabelmat,1)==1; 
                    mask = repmat(stat.negclusterslabelmat(:,maxmaxidx)==negclustidx,1,length(stat.chanlocs));
                else
                    mask = stat.negclusterslabelmat(:,maxmaxidx)==negclustidx;
                end
                topoplot(plotvals,stat.chanlocs, 'maplimits', 'absmax', 'electrodes','off', 'emarker2',{maxchan,'o','green',14,1},...
                    'pmask',mask);
            else
                topoplot(plotvals,stat.chanlocs, 'maplimits', 'absmax', 'electrodes','off', 'emarker2',{maxchan,'o','green',14,1},...
                    'numcontour',0);
            end

            colorbar('FontName',fontname,'FontSize',fontsize);
            title(param.title,'FontName',fontname,'FontSize',fontsize);
        else
            sp=1;
        end

        subplot(sp,1,sp);
        set(gca,'ColorOrder',cat(1,colororder,[0 0 0]));
        hold all;

        if size(stat.diffcond.avg,2)>1
            plot(stat.diffcond.time-stat.timeshift,[stat.diffcond.cond1avg(maxchan,:); stat.diffcond.cond2avg(maxchan,:)]','LineWidth',linewidth*1.5);
        else
            plot(stat.diffcond.time-stat.timeshift,[stat.diffcond.cond1avg, stat.diffcond.cond2avg]','LineWidth',linewidth*1.5);
        end
        %ylim = get(gca,'YLim');
        %ylim = ylim*2;
        ylim = param.ylim;
        set(gca,'YLim',ylim,'XLim',[stat.diffcond.time(1) stat.diffcond.time(end)]-stat.timeshift,'XTick',stat.diffcond.time(1)-stat.timeshift:0.2:stat.diffcond.time(end)-stat.timeshift,...
            'FontName',fontname,'FontSize',fontsize);
        legend(param.legendstrings,'Location',param.legendposition);
        line([stat.diffcond.time(1) stat.diffcond.time(end)]-stat.timeshift,[0 0],'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.60 -0.60],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.45 -0.45],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.30 -0.30],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([-0.15 -0.15],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([    0     0],ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
        line([stat.diffcond.time(maxtime) stat.diffcond.time(maxtime)]-stat.timeshift,ylim,'LineWidth',linewidth,'LineStyle','--','Color','red');
        xlabel('Time (sec) ','FontName',fontname,'FontSize',fontsize);
        ylabel('Amplitude (uV)','FontName',fontname,'FontSize',fontsize);
        box off
        if ~isempty(negclustidx)
            clustwinidx = find(maxval);
            %         rectangle('Position',[stat.time(clustwinidx(1))-stat.timeshift ylim(1) ...
            %             stat.time(clustwinidx(end))-stat.time(clustwinidx(1)) ylim(2)-ylim(1)],'EdgeColor','black','LineStyle','--','LineWidth',linewidth);
            line([stat.time(clustwinidx(1)) stat.time(clustwinidx(end))]-stat.timeshift,[0 0],'Color','blue','LineWidth',8);
            title(sprintf('Cluster @ %.3f sec (t = %.2f, p = %.3f)', stat.diffcond.time(maxtime)-stat.timeshift,stat.negclusters(negclustidx).clusterstat,stat.negclusters(negclustidx).prob),...
                'FontName',fontname,'FontSize',fontsize);
            
            line([stat.time(win(1)) stat.time(win(end))]-stat.timeshift,[0 0],'Color','red','LineWidth',8);
        else
            title(sprintf('%.3f sec', stat.diffcond.time(maxtime)-stat.timeshift),'FontName',fontname,'FontSize',fontsize);
        end
        set(gcf,'Color','white');

    %     export_fig(gcf,[figfile '.eps']);
    end
else
    fprintf('No significant negative clusters found.\n');
end