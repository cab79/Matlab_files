function S = eeg_plot_func(S)
S.func = 'plot';
switch S.plot.type
    case 'topoplot'
        % plot a specific frequency as topoplots
        figure
        datmat = cat(3,S.(S.func).plotdat{:});
        maxdat = max(abs(datmat(:)));
        minmax = [-maxdat,maxdat]+S.plot.scalezero;
        n = 0;
        for c = 1:length(S.(S.func).plotdat)
            dat = S.(S.func).plotdat{c};
            for r = 1:size(dat,2)
                n=n+1;
                subplot(size(dat,2),length(S.(S.func).plotdat),n)
                topoplot(dat(S.(S.func).inclchan,r),S.(S.func).chanlocs(S.(S.func).inclchan),'maplimits',minmax);
                cbar('vert',0,minmax);
                title(S.plot.title{c});
            end
        end
    
%         figure 
%         datmat = cat(3,S.gadata{:}.powspctrm);
%         maxdat = max(abs(datmat(:)));
%         minmax = [-maxdat,maxdat];
%         n = 0;
%         for ii = 1:length(S.gadata);
%             n=n+1;
%             data = S.gadata{ii}.powspctrm;
%             subplot(1,length(S.gadata),n)
%             topoplot(mean(data(:,1),2), chanlocs(S.select.chan{1}),'maplimits',minmax,'electrodes','on'); 
%             cbar('vert',0,minmax);
%         end
%         set(gcf,'Color','white');                                   

    case 'plot_over_freq'
        % plot all frequencies at maximal electrodes
        % find electrodes with maximal power at S.analysistypeed frequency
        dat = cat(3,S.(S.func).plotdat{:});
        [~,maxele] = sort(dat,'descend');
        maxele=maxele(1:S.(S.func).nele);
        %plot
        fplot=[];
        for c = 1:length(S.(S.func).plotdat)
            for fi = 1:length(S.(S.func).select.freqs)
                try
                    fplot(fi,c) = mean(S.(S.func).plotdat{c}(maxele,fi),1);
                catch
                    error('trying to plot more freq than in data - modify S.plot.parts to include ''freqs'' ')
                end
            end
        end
        figure
        plot(S.(S.func).select.freqs,fplot);
        legend(S.plot.legend)
        title(S.plot.title)
        xlabel('Frequency')

    case 'plot_over_time'
        % plot trends for selected freq only
        dat = cat(3,S.(S.func).plotmovdat{:});
        [~,maxele] = sort(dat,'descend');
        maxele=maxele(1:S.(S.func).nele);
        cplot={};
        minlen = inf;
        maxlen = 0;
        for c = 1:length(S.(S.func).plotmovdat)
            cplot{c} = squeeze(mean(mean(S.(S.func).plotmovdat{c}(maxele,:,:),1),2));
            minlen = min(minlen,length(cplot{c}));
            maxlen = max(maxlen,length(cplot{c}));
        end
        if 0
            % equate lengths
            for c = 1:length(S.(S.func).plotmovdat)
                cplot{c} = cplot{c}(1:minlen);
            end
        else % pad
            for c = 1:length(S.(S.func).plotmovdat)
                if length(cplot{c})<maxlen
                    cplot{c}(end+1:maxlen) = nan;
                end
            end
        end
        fplot = cat(2,cplot{:});
        figure
        plot(1:size(fplot,1),fplot);
        legend(S.plot.legend)
        title(S.plot.title)
        xlabel(['Time in trials)'])

    case 'timtopo'
        for ga = 1:length(S.(S.func).gadata);
            plotdata = S.(S.func).gadata{ga}.avg;
            [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
            [~, maxmaxidx] = max(maxval);
            plottime = S.(S.func).gadata{ga}.time(maxidx(maxmaxidx));
            if plottime == S.(S.func).gadata{ga}.time(end)
                plottime = S.(S.func).gadata{ga}.time(end-1);
            end
            figure('Name','grand average');
            timtopo(plotdata,S.(S.func).chanlocs(S.(S.func).select.chan{1}),...
                'limits',[S.(S.func).gadata{ga}.time(1) S.(S.func).gadata{ga}.time(end)]*1000,...
                'plottimes',plottime*1000);
            set(gcf,'Color','white');
        end
        
    case 'plottopo'
        figure
        plottopo([S.(S.func).gadata{ga}.avg,S.(S.func).gadata{ga}.avg-S.(S.func).gadata{ga}.var,S.(S.func).gadata{ga}.avg+S.(S.func).gadata{ga}.var], 'chanlocs', S.(S.func).chanlocs(S.(S.func).select.chan{1}), 'frames', length(S.(S.func).gadata{ga}.time), ...
          'limits', [S.(S.func).gadata{ga}.time(1) S.(S.func).gadata{ga}.time(end) 0 0]*1000, 'title', '', 'colors', { 'k' 'r--' 'r--' }, ...
          'chans', 0, 'legend', {'avg','std','std'}, 'regions', {}, 'ylim', []);
        set(gcf,'Color','white');
    
end
    
