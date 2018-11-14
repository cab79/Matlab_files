function S=gfp_analysis(S)

dbstop if error
S.func = 'gfp';
if isfield(S.(S.func),'GFP')
    S.(S.func) = rmfield(S.(S.func),'GFP');
end

% select channels
S=select_chans(S);


% if there is no data within S already (e.g. from just running FT_grand_average), then load it
%if ~isfield(S,'gadata')
    switch S.(S.func).select.datatype
        case 'ERP'
            S.path.file = S.path.erp;
        case 'TF'
            S.path.file = S.path.freq;
        case 'Freq'
            error('GFP analysis not possible with Freq data, only ERP and TF data')
    end
    
    files = dir(fullfile(S.path.file,'*grandavg*mat'));
    for ga = 1:length(files)
        load(fullfile(S.path.file,files(ga).name));
        try
            S.(S.func).gadata{ga} = gadata.gavg;
        catch
            S.(S.func).gadata{ga} = gadata;
        end
        [~,S.(S.func).ganame{ga},~] = fileparts(files(ga).name);
    end
%end

% normalise by common average and GFP mean over time
for ga = 1:length(S.(S.func).gadata)
    if S.(S.func).GAnorm.commonavg
        for e = 1:size(S.(S.func).gadata{ga}.avg,1)
            S.(S.func).gadata{ga}.avg(e,:) = S.(S.func).gadata{ga}.avg(e,:) - nanmean(S.(S.func).gadata{ga}.avg,1);
        end
    end
    if S.(S.func).GAnorm.gfp
        S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.avg/nanmean(std(S.(S.func).gadata{ga}.avg,1));
    end
    if S.(S.func).GArmbase
        correcttimes = dsearchn(S.(S.func).gadata{ga}.time',S.(S.func).epoch.basewin');
        basemean = mean(S.(S.func).gadata{ga}.avg(:,correcttimes(1):correcttimes(2)),2);
        S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.avg - repmat(basemean,1,size(S.(S.func).gadata{ga}.avg,2));
    end
end

% smooth gavg for peak identification
for ga = 1:length(S.(S.func).gadata)
    for e = 1:size(S.(S.func).gadata{ga}.avg,1)
        S.(S.func).gadata{ga}.avg(e,:) = smooth(squeeze(S.(S.func).gadata{ga}.avg(e,:)),30,'lowess');
    end
end

% identify time bins
gatimes = S.(S.func).gadata{1}.time;
basetimesGFP = intersect(find(gatimes<S.(S.func).epoch.basewin(2)), find(gatimes>=S.(S.func).epoch.basewin(1)));
stimtimesGFP = find(gatimes>S.(S.func).epoch.basewin(2));
timesRMSE = gatimes(2:end);
basetimesRMSE = intersect(find(timesRMSE<S.(S.func).epoch.basewin(2)), find(timesRMSE>=S.(S.func).epoch.basewin(1)));
stimtimesRMSE = find(timesRMSE>S.(S.func).epoch.basewin(2));

% GFP and CIs across subjects
for ga = 1:length(S.(S.func).gadata)
    GFP = std(S.(S.func).gadata{ga}.avg,1);
    GFP_CI = repmat(mean(GFP(basetimesGFP))+std(GFP(basetimesGFP)),length(GFP),1);
    GFP_base = repmat(mean(GFP(basetimesGFP)),length(GFP),1);
    figure
    yyaxis left
    plot(gatimes,GFP);
    hold on
    plot(gatimes,GFP_CI,'b');
    hold on
    plot(gatimes,GFP_base,'b--');

    % GFP peaks
    [peaks_GFP peaklocs_GFP] = findpeaks_Rbase(GFP,'MinPeakHeight',unique(GFP_CI),'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    peaklocgatimes_GFP = gatimes(peaklocs_GFP);
    hold on;
    plot(peaklocgatimes_GFP,peaks_GFP,'v','MarkerFaceColor','b');
    % GFP troughs
    [tro_GFP trolocs_GFP] = findpeaks_Lbase(-GFP,'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    trolocgatimes_GFP = gatimes(trolocs_GFP);
    tro_GFP = -tro_GFP;
    hold on;
    plot(trolocgatimes_GFP,tro_GFP,'s','MarkerFaceColor','b');

    % RMSE
    RMSE = [];
    for t = 2:size(S.(S.func).gadata{ga}.avg,2)
        y = S.(S.func).gadata{ga}.avg(:,t);
        yhat = S.(S.func).gadata{ga}.avg(:,t-1);
        RMSE(t-1) = sqrt(mean((y - yhat).^2));
    end
    RMSE_CI = repmat(mean(RMSE(basetimesRMSE))+std(RMSE(basetimesRMSE)),length(RMSE),1);
    RMSE_base = repmat(mean(RMSE(basetimesRMSE)),length(RMSE),1);
    hold on
    yyaxis right
    plot(gatimes(2:end),RMSE,'r');
    hold on
    plot(gatimes(2:end),RMSE_CI,'r');
    hold on
    plot(gatimes(2:end),RMSE_base,'r--');

    % RMSE peaks
    [peaks_RMSE peaklocs_RMSE] = findpeaks_Rbase(RMSE,'MinPeakHeight',unique(RMSE_CI),'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    peaklocs_RMSE = peaklocs_RMSE+1;
    hold on;
    plot(gatimes(peaklocs_RMSE),peaks_RMSE,'v','MarkerFaceColor','r');
    % RMSE troughs
    [tro_RMSE trolocs_RMSE] = findpeaks_Lbase(-RMSE,'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    trolocs_RMSE = trolocs_RMSE+1;
    tro_RMSE = -tro_RMSE;
    hold on;
    plot(gatimes(trolocs_RMSE),tro_RMSE,'s','MarkerFaceColor','r');

    % find peak limits
    peaklocs_RMSE=peaklocs_RMSE(peaklocs_RMSE>stimtimesRMSE(1));
    limit1=[];
    limit2=[];
    if ~isempty(peaklocs_RMSE)
        for i = 1:length(peaklocs_RMSE)
            limit1=[limit1 peaklocs_RMSE(i)+1];
            if i<length(peaklocs_RMSE)
                next_tros = trolocs_RMSE(intersect(find(trolocs_RMSE>peaklocs_RMSE(i)),find(trolocs_RMSE<peaklocs_RMSE(i+1))));
            else
                next_tros = trolocs_RMSE(find(trolocs_RMSE>peaklocs_RMSE(i)));
            end
            if ~isempty(next_tros)
                limit2=[limit2 next_tros(end)];
            else
                % NOT SURE WHICH IS BEST HERE:
                %limit1 = limit1(1:end-1);
                %continue
                % OR:
                limit2 = stimtimesRMSE(end);
            end
            S.(S.func).GFP(ga).limits_all{i}=[limit1(end):limit2(end)];
        end
        S.(S.func).GFP(ga).limits_all=S.(S.func).GFP(ga).limits_all(~cellfun('isempty',S.(S.func).GFP(ga).limits_all));
    else
        S.(S.func).GFP(ga).limits_all = {};
    end

    if isempty(S.(S.func).GFP(ga).limits_all)
        no_peaks_temp=0;
    else
        no_peaks_temp=length(S.(S.func).GFP(ga).limits_all);
    end

    %TOPOSIM
    if no_peaks_temp>1
        TOPO = [];
        for n = 1:no_peaks_temp
            TOPO(:,n) = mean(S.(S.func).gadata{ga}.avg(:,S.(S.func).GFP(ga).limits_all{n}),2);
            %subplot(length(subjects),no_peaks,(s-1)*no_peaks+n); topoplot(TOPO(:,s,n), EEG.chanlocs);
        end
        TOPOSIM = pdist(TOPO','cosine');
        TOPOSIM=squareform(TOPOSIM);

        for thresh=0;
            test=TOPOSIM.*(TOPOSIM<thresh);
            test = tril(test);
            for t = 1:size(test,2)
                if t<size(test,2)-1
                    test(t+2:end,t)=0; % only consider adjacent connections
                end
            end
            test = test'+test;
            bins = conncomp(graph(test));
            test2 = bins(2:end)-bins(1:end-1);
        end

        no_peaks=max(bins);
        for n = 1:no_peaks
            bin_idx = find(bins==n);
            no_temp_peaks = length(bin_idx);
            limits_temp_np = S.(S.func).GFP(ga).limits_all(bin_idx);
            S.(S.func).GFP(ga).limits{n}=[];
            for tp = 1:no_temp_peaks
                S.(S.func).GFP(ga).limits{n}= [S.(S.func).GFP(ga).limits{n} limits_temp_np{tp}];
            end
        end
        figure
        for n = 1:no_peaks
            subplot(1,no_peaks,n); topoplot(mean(S.(S.func).gadata{ga}.avg(:,S.(S.func).GFP(ga).limits{n}),2), S.(S.func).chanlocs(S.(S.func).inclchan));
        end

        % plot grand average with time bins
        plotdata = S.(S.func).gadata{ga}.avg;
        figure('Name','grand average');
        plot(gatimes,plotdata');
        set(gcf,'Color','white');
        cm=colormap(jet(no_peaks));
        hold on
        for p = 1:no_peaks
            binidx = find(bins==p);
            for n = 1:length(binidx)
                line([gatimes(S.(S.func).GFP(ga).limits_all{binidx(n)}(1)) gatimes(S.(S.func).GFP(ga).limits_all{binidx(n)}(end))],[0 0],'color',cm(p,:),'LineWidth',8);
            end
        end
        
    else
        no_peaks=0;
        bins=0;
    end
    S.(S.func).GFP(ga).bins_all=bins;
end

% identify GFPs for each subject and condition
if ~isfield(S.ga,'gafiles')
    error('run eeg_grand_average before this function to ensure S.(S.func).gafiles is present');
end
for ga = 1:length(S.(S.func).gadata)
    
    % create tables. first column is filename
    T_GFP=table;
    T_POW=table;
    T_LAT=table;
    clear C_GFP C_POW C_LAT
    C_GFP = ['Filename'; S.ga.gafiles{ga}'];
    C_POW=C_GFP;
    C_LAT=C_GFP;
    
    for f = 1:length(S.ga.gafiles{ga}) 
        
        % load data for each file
        load(fullfile(S.path.file,S.ga.gafiles{ga}{f}));
        
        switch S.(S.func).select.datatype
            case 'ERP'
                % select events
                if isstruct(tldata)
                    data_all = {tldata};
                elseif iscell(tldata)
                    data_all = tldata(S.(S.func).select.events);
                end
            case 'TF'
                % select events
                if isstruct(fdata)
                    data_all = {fdata};
                elseif iscell(fdata)
                    data_all = fdata(S.(S.func).select.events);
                end
                % select frequency using S.(S.func).fn 
                % create .trial and .avg
        end
        
        %for each event type
        col=1;
        for ev=1:length(data_all)
            dat = data_all{ev};
        
            % for each peak, extract data
            no_peaks = length(S.(S.func).GFP(ga).limits);
            for p = 1:no_peaks
                
                GFP = std(dat.avg(S.(S.func).inclchan,S.(S.func).GFP(ga).limits{p}),[],1);
                [~, GFPmaxLOC] = max(GFP);
                GFP = sort(GFP(:), 'descend');
                switch S.(S.func).GFPnorm.type
                    case 'GFPnorm'
                        GFP=GFP/mean(squeeze(std(dat.avg,[],1)));
                    case 'GFPpriornorm1'
                        GFP=GFP/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{p}(1)-1),[],1)));
                    case 'GFPpriornorm2'
                        if p==1; 
                            GFP=GFP/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{p}(1)-1),[],1)));
                        else
                            GFP=GFP/mean(squeeze(std(dat.avg(:,S.(S.func).GFP(ga).limits{p-1}),[],1)));
                        end
                    case 'GFPbase'
                        GFP=GFP-mean(squeeze(std(dat.avg(:,basetimesGFP),[],1)));
                    case 'GFPp1norm'
                        GFP=GFP/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{1}),[],1)));
                    case 'lognorm'
                        GFP=log(GFP);
                end

                POW = mean(dat.avg(S.(S.func).inclchan,S.(S.func).GFP(ga).limits{p}),2);
                POW = sort(POW(:), 'descend');
                %if posnegpeak(p)==1
                %    POW = sort(POW(:), 'descend');
                %elseif posnegpeak(p)==2
                %    POW = sort(POW(:), 'ascend');
                %end

                switch S.(S.func).GFPnorm.type
                    case 'GFPnorm'
                        POW=POW/mean(squeeze(std(dat.avg,[],1)));
                    case 'GFPpriornorm1'
                        POW=POW/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{p}(1)-1),[],1)));
                    case 'GFPpriornorm2'
                        if p==1; 
                            POW=POW/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{p}(1)-1),[],1)));
                        else
                            POW=POW/mean(squeeze(std(dat.avg(:,S.(S.func).GFP(ga).limits{p-1}),[],1)));
                        end
                    case 'GFPbase'
                        POW=POW-mean(squeeze(std(dat.avg(:,basetimesGFP),[],1)));
                    case 'GFPp1norm'
                        POW=POW/mean(squeeze(std(dat.avg(:,1:S.(S.func).GFP(ga).limits{1}),[],1)));
                    case 'lognorm'
                        POW=log(POW);
                end

                GFP = mean(GFP(1:ceil(length(GFP)*S.(S.func).topXpercent/100)));
                POW = mean(POW(1:ceil(length(GFP)*S.(S.func).topXpercent/100)));
                lat_range = timesRMSE(S.(S.func).GFP(ga).limits{p});
                LAT = lat_range(GFPmaxLOC)*1000;
                
                col = col+1;
                % add data to table
                header = ['event' num2str(ev) '_peak' num2str(p)];
                C_GFP{1,col} = header;
                C_POW{1,col} = header;
                C_LAT{1,col} = header;
                C_GFP{f+1,col} = GFP;
                C_POW{f+1,col} = POW;
                C_LAT{f+1,col} = LAT;

            end
        end
    end
    T_GFP=cell2table(C_GFP(2:end,:));
    T_GFP.Properties.VariableNames = C_GFP(1,:);
    T_POW=cell2table(C_POW(2:end,:));
    T_POW.Properties.VariableNames = C_POW(1,:);
    T_LAT=cell2table(C_LAT(2:end,:));
    T_LAT.Properties.VariableNames = C_LAT(1,:);
    datename = datestr(now,30); 
    cd(S.path.file)
    writetable(T_GFP,['gfp_results_' S.(S.func).ganame{ga} '_' datename '.xlsx']);
    writetable(T_POW,['pow_results_' S.(S.func).ganame{ga} '_' datename '.xlsx']);
    writetable(T_LAT,['lat_results_' S.(S.func).ganame{ga} '_' datename '.xlsx']);
end
