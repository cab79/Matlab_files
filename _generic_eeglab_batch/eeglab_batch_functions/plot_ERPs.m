function S=plot_ERPs(S)

dbstop if error

S.func = 'ploterp';
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end

% GET FILE LIST
S = getfilelist(S);

% select channels
S=select_chans(S);

% show outlier-ness
if exist(fullfile(S.path.file,'Outliers.mat'),'file');
    load(fullfile(S.path.file,'Outliers.mat'));
    subs = S.(S.func).designmat(2:end,find(strcmp(S.(S.func).designmat(1,:),'subjects')));
end

for i = 1:length(S.(S.func).filelist)
    
    file = S.(S.func).filelist{i};
    load(file)
    noemp = find(~cellfun(@isempty,tldata));
    tldata = tldata(noemp);
    avgdata = tldata{1};
    datstruct = cell2mat(tldata);
    datmat = cat(3,datstruct(:).avg);
    avgdata.avg = mean(datmat,3);
    
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    cfg = [];
    cfg.layout = S.ploterp.layout;
    cfg.ylim = [prctile(datmat(:),0.1),prctile(datmat(:),99.9)];%S.ploterp.ylim;
    for d = 1:length(tldata)
        cfg.dataname{d} = ['num trials: ' num2str(size(tldata{d}.trial,1))];
    end
    ft_multiplotER_cab(cfg, tldata);
    title(file)
    
    f2 = figure('units','normalized','outerposition',[0 0.8 0.1*length(S.ploterp.times) 0.2]);
    for t = 1:length(S.ploterp.times)
        subplot(1,length(S.ploterp.times),t);
        tim = dsearchn(avgdata.time',S.ploterp.times{t}');
        topoplot(mean(avgdata.avg(:,tim(1):tim(2)),2),S.(S.func).chanlocs(S.(S.func).inclchan));
    end
    %cfg = [];                            
    %cfg.xlim = [0.3 0.5];                
    %cfg.zlim = [0 6e-14];                
    %cfg.layout = S.ploterp.layout;            
    %cfg.parameter = 'avg'; % the default 'avg' is not present in the data
    %figure; ft_topoplotER(cfg,avgdata); 
    
    % show outlier-ness
    if exist('outlist','var');
        max([outlist{:,2}])
        %outval = round(100*outlist{strcmp(outlist(:,1),subs{i}),2}/max([outlist{:,2}]));
        outval = round(invprctile([outlist{:,2}],outlist{strcmp(outlist(:,1),subs{i}),2}));
    end
    
    % good, bad?
    S.(S.func).qual{i,1} = file;
    S.(S.func).qual{i,2} = MFquestdlg([0.9 0.3],['outlier percentile: ' num2str(outval)],'Data quality','Good','So-so','Bad','So-so')
    close(f); close(f2)
    qual=S.(S.func).qual;
    save(fullfile(S.path.file,'data_quality'),'qual')
    
    if strcmp(S.(S.func).qual{i,2},'')
        return
    end

end

