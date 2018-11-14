function S=plot_ERPs(S,varargin)

dbstop if error

if nargin>1
    type = varargin{1};
else
    type = 'subject';
end

S.func = 'ploterp';
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end

% select channels
S=select_chans(S);

switch type
    case 'subject'

    % GET FILE LIST
    S = getfilelist(S);

    % show outlier-ness
    if exist(fullfile(S.path.file,'Outliers.mat'),'file')
        load(fullfile(S.path.file,'Outliers.mat'));
        subs = S.(S.func).designmat(2:end,find(strcmp(S.(S.func).designmat(1,:),'subjects')));
    end

    for i = 1:length(S.(S.func).filelist)

        file = S.(S.func).filelist{i};
        load(file)
        
        % select events
        if iscell(S.ploterp.select.events)
            for ev = 1:length(S.ploterp.select.events)
                newdata{ev} = tldata{1};
                selectdata = tldata(S.ploterp.select.events{ev});
                datstruct = cell2mat(selectdata(~cellfun(@isempty,selectdata)));
                datmat = cat(1,datstruct(:).trial);
                newdata{ev}.trial = datmat;
                newdata{ev}.avg = squeeze(mean(newdata{ev}.trial,1));
            end
            tldata=newdata;
        else
            tldata = tldata(S.ploterp.select.events);
        end
        
        noemp = find(~cellfun(@isempty,tldata));
        tldata = tldata(noemp);
        avgdata = tldata{1};
        datstruct = cell2mat(tldata);
        datmat = cat(3,datstruct(:).avg);
        avgdata.avg = mean(datmat,3);

        [f,f2]=plotmulti(S,datmat,tldata,file,avgdata)

        % show outlier-ness
        if exist('outlist','var');
            max([outlist{:,2}])
            %outval = round(100*outlist{strcmp(outlist(:,1),subs{i}),2}/max([outlist{:,2}]));
            outval = round(invprctile([outlist{:,2}],outlist{strcmp(outlist(:,1),subs{i}),2}));
        else
            outval = 0;
        end

        % good, bad?
        S.(S.func).qual{i,1} = file;
        S.(S.func).qual{i,2} = MFquestdlg([0.9 0.3],['outlier percentile: ' num2str(outval)],'Data quality','Good','So-so','Bad','So-so')
        if isempty(S.(S.func).qual{i,2})
            return
        end
        close(f); close(f2)
        qual=S.(S.func).qual;
        save(fullfile(S.path.file,'data_quality'),'qual')

        if strcmp(S.(S.func).qual{i,2},'')
            return
        end

    end
    
    case 'grandavg'
        
        if ~isfield(S.ga,'ganame')
            S.ga.ganame = {'grandavg.mat'};
        end
        
        for ga = 1:length(S.ga.ganame)
            
            load(fullfile(S.path.file,S.ga.ganame{ga}));
            avgdata = gadata.gavg;
            tldata = gadata.events;
            datstruct = cell2mat(tldata);
            datmat = cat(3,datstruct(:).avg);

            title = 'grand average';
            [f,f2]=plotmulti(S,datmat,tldata,title,avgdata)
            
            [f3,f4]=plotvar(S,datmat,tldata,title,avgdata)

        end
        
end

function [f,f2] = plotmulti(S,datmat,tldata,file,avgdata)
f=figure('units','normalized','outerposition',[0 0 1 1]);
cfg = [];
cfg.layout = S.ploterp.layout;
cfg.ylim = [prctile(datmat(:),0.1),prctile(datmat(:),99.9)];%S.ploterp.ylim;
if ~isfield(S.ploterp,'event_labels') || isempty(S.ploterp.event_labels)
    temp=1:length(tldata);
    labels = cellstr(num2str(temp(:)))';
else
    labels = S.ploterp.event_labels;
end
if isfield(tldata{1},'trial')
    for d = 1:length(tldata)
        cfg.dataname{d} = ['event: ' labels{d} ', num trials: ' num2str(size(tldata{d}.trial,1))];
    end
else
    for d = 1:length(tldata)
        cfg.dataname{d} = ['event: ' labels{d}];
    end
end
ft_multiplotER_cab(cfg, tldata);
title(file)

f2 = figure('units','normalized','outerposition',[0 0.7 0.1*length(S.ploterp.times) 0.3]);
for t = 1:length(S.ploterp.times)
    subplot(1,length(S.ploterp.times),t);
    tim = dsearchn(avgdata.time',S.ploterp.times{t}');
    topoplot(mean(avgdata.avg(:,tim(1):tim(2)),2),S.(S.func).chanlocs(S.(S.func).inclchan));
    title([num2str(S.ploterp.times{t}(1)) ':' num2str(S.ploterp.times{t}(2))])
end
%cfg = [];                            
%cfg.xlim = [0.3 0.5];                
%cfg.zlim = [0 6e-14];                
%cfg.layout = S.ploterp.layout;            
%cfg.parameter = 'avg'; % the default 'avg' is not present in the data
%figure; ft_topoplotER(cfg,avgdata); 


function [f3,f4] = plotvar(S,datmat,tldata,file,avgdata)

% settings
lat_ms = [268]; %ms
pos_neg_peak = [1,-1]; %1 or -1
num_ele = 3;
plotidx=[1 2];
col=colormap(jet(length(tldata)));

for i = 1:length(S.ploterp.times)
    tim{i} = dsearchn(avgdata.time',S.ploterp.times{i}');
end

for p = 1:length(pos_neg_peak)

    for i = 1:length(S.ploterp.times)
        
        % identify peak electrodes
        if pos_neg_peak(p)==1
            [~,tps] = sort(mean(avgdata.avg(:,tim{i}(1):tim{i}(2)),2),'descend');
            tp{i} = tps(1:num_ele);
        elseif pos_neg_peak(p)==-1
            [~,tps] = sort(mean(avgdata.avg(:,tim{i}(1):tim{i}(2)),2),'ascend');
            tp{i} = tps(1:num_ele);
        end
        
        f3(p,i)=figure('units','normalized','outerposition',[0 0 1 1]);
        ts=avgdata.time;
        ptimes = 1:length(avgdata.time);
        hold all
        st=avgdata.time(tim{i}(1):tim{i}(2));
        fill([st, fliplr(st)], [ones(1,length(st))*100, fliplr(ones(1,length(st))*-100)], ...
        [0.5 0.5 0.5], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
        %plot([lat_ms lat_ms],[-10 10],'k');
        lowest=Inf;
        highest=-Inf;
        for ii = 1:length(tldata)
            peakdata = tldata{ii}.avg;
            vardata = tldata{ii}.var;
            ERP = mean(peakdata(tp{i},:),1);
            VAR = mean(sqrt(vardata(tp{i},:)),1);
            %nsub = length(tldata(ii,:));
            %SEM = VAR/sqrt(nsub);               % Standard Error
            %tscore = -tinv(0.025,nsub-1);      % T-Score
            %CI = tscore*SEM;                      % Confidence Intervals
            upper = ERP(ptimes)+VAR(ptimes);
            lower = ERP(ptimes)-VAR(ptimes);
            lowest = min(min(lower),lowest);
            highest = max(max(upper),highest);
            fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
            col(ii,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
            plot(ts,ERP(ptimes),'color',col(ii,:)); 
        end
        ylabel('Amplitude, uV') % label for y axis
        xlabel('Time (ms)') % label for x axis
        set(gca,'FontSize',15);
        set(findall(gcf,'type','text'),'FontSize',15);
        ylim([lowest highest]);
        hold off
    end

    f4(p) = figure('units','normalized','outerposition',[0 0 0.1*length(tldata) 0.2*length(S.ploterp.times)]);
    plotchans=1:length(S.(S.func).chanlocs);
    iii=0;
    for i = 1:length(S.ploterp.times)
        for ii = 1:length(tldata);
            iii = iii+1;
            peakdata = tldata{ii}.avg;
            [~,markchans] = intersect(plotchans,tp{i});
            subplot(length(tldata),length(S.ploterp.times),iii); 
            topoplot(mean(peakdata(:,tim{i}(1):tim{i}(2)),2), S.(S.func).chanlocs(S.(S.func).inclchan),'maplimits',[-4 4],'electrodes','on','plotchans',plotchans,'emarker2',{markchans,'o','w',3,1}); colorbar
        end
    end
end
