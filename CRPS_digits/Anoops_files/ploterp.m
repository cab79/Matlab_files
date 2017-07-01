function ploterp(statmode,subjinfo,condlist,varargin)

loadpaths

timeshift = 0; %milliseconds

param = finputcheck(varargin, { 'ylim', 'real', [], [-5 5]; ...
    'subcond', 'string', {'on','off'}, 'off'; ...
    'topowin', 'real', [], []; ...
    'legendstrings', 'cell', {}, condlist; ...
    });

%% SELECTION OF SUBJECTS AND LOADING OF DATA

loadsubj;

if strcmp(statmode,'trial')
    %%%% perform single-trial statistics
    subjlist = {subjinfo};
    subjcond = condlist;
    
elseif strcmp(statmode,'cond')
    %%%% perform within-subject statistics
    subjlist1 = subjlists{subjinfo(1)};
    subjlist2 = subjlists{subjinfo(2)};
    subjlist = cat(2,subjlist1,subjlist2);
    subjcond = repmat(condlist,length(subjlist),1);

elseif strcmp(statmode,'subj')
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
end

numsubj = length(subjlist);
numcond = size(subjcond,2);

conddata = cell(numsubj,numcond);

%% load and prepare individual subject datasets

for s = 1:numsubj
    
    % %     rereference
    %     EEG = rereference(EEG,1);
    
        %%%%% baseline correction relative to 5th tone
%         bcwin = [-200 0];
%         bcwin = bcwin+timeshift;
%         EEG = pop_rmbase(EEG,bcwin);
        %%%%%
    
    
    for c = 1:numcond
        EEG = pop_loadset('filename', sprintf('%s.set', subjlist{s,c}), 'filepath', filepath);
        %EEG = sortchan(EEG);
        % THIS ASSUMES THAT ALL DATASETS HAVE SAME NUMBER OF ELECTRODES
        if s == 1
            chanlocs = EEG.chanlocs;
            erpdata = zeros(EEG.nbchan,EEG.pnts,numcond,numsubj);
        end
        
        conddata{s,c} = EEG;
        
%         saveEEG = conddata{s,c};
%         saveEEG.data = mean(conddata{s,c}.data,3);
%         saveEEG.setname = sprintf('%s_%s_%s',statmode,subjlist{s},subjcond{s,c});
%         saveEEG.filename = [saveEEG.setname '.set'];
%         saveEEG.trials = 1;
%         saveEEG.event = saveEEG.event(1);
%         saveEEG.event(1).type = saveEEG.setname;
%         saveEEG.epoch = saveEEG.epoch(1);
%         saveEEG.epoch(1).eventtype = saveEEG.setname;
%         pop_saveset(saveEEG,'filepath',filepath,'filename',saveEEG.filename);
        
%         if (strcmp(statmode,'trial') || strcmp(statmode,'cond')) && c == numcond
%             if conddata{s,1}.trials > conddata{s,2}.trials
%                 fprintf('Equalising trials in condition %s.\n',subjcond{s,1});
%                 randtrials = randperm(conddata{s,1}.trials);
%                 conddata{s,1} = pop_select(conddata{s,1},'trial',randtrials(1:conddata{s,2}.trials));
%             elseif conddata{s,2}.trials > conddata{s,1}.trials
%                 fprintf('Equalising trials in condition %s.\n',subjcond{s,2});
%                 randtrials = randperm(conddata{s,2}.trials);
%                 conddata{s,2} = pop_select(conddata{s,2},'trial',randtrials(1:conddata{s,1}.trials));
%             end
%         end
        
        erpdata(:,:,c,s) = mean(conddata{s,c}.data,3);
    end
end

if strcmp(statmode,'subj')
    if numcond == 2
        erpdata = erpdata(:,:,1,:) - erpdata(:,:,2,:);
        condlist = {sprintf('%s-%s',condlist{1},condlist{3}),sprintf('%s-%s',condlist{2},condlist{4})};
    end
    if strcmp(param.subcond,'on')
        erpdata = mean(erpdata(:,:,:,1:numsubj1),4) - mean(erpdata(:,:,:,numsubj1+1:end),4);
        condlist = {sprintf('%s-%s',condlist{1},condlist{2})};
    else
        erpdata = cat(3, mean(erpdata(:,:,:,1:numsubj1),4), mean(erpdata(:,:,:,numsubj1+1:end),4));
    end
else
    if strcmp(statmode,'cond') && numcond == 3
        erpdata(:,:,1,:) = erpdata(:,:,1,:) - erpdata(:,:,3,:);
        erpdata(:,:,2,:) = erpdata(:,:,2,:) - erpdata(:,:,3,:);
        erpdata = erpdata(:,:,[1 2],:);
        numcond = 2;
        condlist = {sprintf('%s-%s',condlist{1},condlist{3}),sprintf('%s-%s',condlist{2},condlist{3})};
    end
    
    if numcond == 2 && strcmp(param.subcond,'on')
        for s = 1:size(erpdata,4)
            erpdata(:,:,3,s) = erpdata(:,:,1,s) - erpdata(:,:,2,s);
        end
        condlist = cat(2,condlist,{sprintf('%s-%s',condlist{1},condlist{2})});
        
        erpdata = erpdata(:,:,3,:);
        condlist = condlist(3);
    end
    erpdata = mean(erpdata,4);
end

%% PLOTTING

for c = 1:size(erpdata,3)
    plotdata = erpdata(:,:,c);
    
    if isempty(param.topowin)
        param.topowin = [0 EEG.times(end)-timeshift];
    end
    latpnt = find(EEG.times-timeshift >= param.topowin(1) & EEG.times-timeshift <= param.topowin(2));
    [maxval, maxidx] = max(abs(plotdata(:,latpnt)),[],2);
    [~, maxmaxidx] = max(maxval);
    plottime = EEG.times(latpnt(1)-1+maxidx(maxmaxidx));
    if plottime == EEG.times(end)
        plottime = EEG.times(end-1);
    end
    
    %plot ERP data
    figure('Name',condlist{c});
    timtopo(plotdata,chanlocs,...
        'limits',[EEG.times(1)-timeshift EEG.times(end)-timeshift, param.ylim],...
        'plottimes',plottime-timeshift);
    set(gcf,'Color','white');
    
    %saveEEG = pop_select(EEG,'point',[latpnt(1) latpnt(end)+1]);
    %saveEEG.data = plotdata(:,latpnt);
%     saveEEG = EEG;
%     saveEEG.data = plotdata;
%     if strcmp(statmode,'cond') || strcmp(statmode,'trial')
%         saveEEG.setname = sprintf('%s_%s_%s',statmode,num2str(subjinfo),condlist{c});
%     elseif strcmp(statmode,'subj')
%         saveEEG.setname = sprintf('%s_%s_%s',statmode,num2str(subjinfo(c)),condlist{c});
%     end
%     saveEEG.filename = [saveEEG.setname '.set'];
%     saveEEG.trials = 1;
%     saveEEG.event = saveEEG.event(1);
%     saveEEG.event(1).type = saveEEG.setname;
%     saveEEG.epoch = saveEEG.epoch(1);
%     saveEEG.epoch(1).eventtype = saveEEG.setname;
%     pop_saveset(saveEEG,'filepath',filepath,'filename',saveEEG.filename);
end
% 
% fontname = 'Helvetica';
% fontsize = 20;
% linewidth = 2;
% plotchan = 'E10';
% plotchan = find(strcmp(plotchan,{EEG.chanlocs.labels}));
% 
% colorlist = {
%     'local standard'    [0         0    1.0000]
%     'local deviant'     [0    0.5000         0]
%     'global standard'   [1.0000    0         0]
%     'global deviant'    [0    0.7500    0.7500]
%     'inter-aural dev.'  [0.7500    0    0.7500]
%     'inter-aural ctrl.' [0.7500    0.7500    0]
%     'attend tones'      [0    0.5000    0.5000]
%     'attend sequences'  [0.5000    0    0.5000]
%     'attend visual'     [0    0.2500    0.7500]
%     'early glo. std.'   [0.5000    0.5000    0]
%     'late glo. std.'    [0.2500    0.5000    0]
%     };
% 
% curcolororder = get(gca,'ColorOrder');
% colororder = zeros(length(param.legendstrings),3);
% for str = 1:length(param.legendstrings)
%     cidx = strcmp(param.legendstrings{str},colorlist(:,1));
%     if sum(cidx) == 1
%         colororder(str,:) = colorlist{cidx,2};
%     else
%         colororder(str,:) = curcolororder(str,:);
%     end
% end
% 
% newfig = false;
% 
% if newfig
% figure;
% figpos = get(gcf,'Position');
% set(gcf,'Position',[figpos(1) figpos(2) figpos(3)*2 figpos(4)]);
% end
% 
% set(gca,'ColorOrder',colororder);
% hold all
% plot(EEG.times-timeshift,squeeze(erpdata(plotchan,:,:)),'LineWidth',linewidth*1.5);
% for c = 1:size(erpdata,3)
%         condfit = polyfit(EEG.times(latpnt),erpdata(plotchan,latpnt,c),1);
%         plot(EEG.times-timeshift,polyval(condfit,EEG.times),'LineWidth',linewidth,'LineStyle','--','Color',colororder(c,:));
% end
% 
% if newfig
%     set(gca,'XLim',[EEG.times(1) EEG.times(end)]-timeshift,'YLim',param.ylim,'XTick',EEG.times(1)-timeshift:200:EEG.times(end)-timeshift,...
%         'FontName',fontname,'FontSize',fontsize);
%     line([EEG.times(1) EEG.times(end)]-timeshift,[0 0],'LineWidth',linewidth,'Color','black','LineStyle',':');
%     line([0 0],param.ylim,'LineWidth',linewidth,'Color','black','LineStyle',':');
%     xlabel('Time relative to 5th tone (sec) ','FontName',fontname,'FontSize',fontsize);
%     ylabel('Amplitude (uV)','FontName',fontname,'FontSize',fontsize);
%     box on
%     legend(param.legendstrings,'Location','NorthWest');
% else
%     [~,~,~,curlegend] = legend;
%     legend(cat(2,curlegend,param.legendstrings),'Location','NorthWest');
% end
% set(gcf,'Color','white');