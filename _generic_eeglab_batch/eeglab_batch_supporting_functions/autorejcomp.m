function [EEG temp detrendcomps] = autorejcomp(EEG,option,varargin)

%varargin(1) = [start,end] time window in s to search for peaks/troughs
%varargin(2) = threshold % of trialwise consistency for accepting component
%OR ratio of peaks vs baseline on variance measures
%varargin(3) = maxmimum component to include
%varargin(4) = baseline

%EEG.reject.gcompreject = ones(1,length(EEG.icachansind));

        temp=[];
if option==1
    EEG.reject.gcompreject = ones(1,size(EEG.icawinv,2));
elseif option==2 && ~isempty(varargin)
    lats = varargin{1}*1000;
    for i = 1:size(EEG.icawinv,2)
        icaacttmp  = squeeze(eeg_getdatact(EEG, 'component', i)); % provides dp x trials
        offset     = nan_mean(icaacttmp(EEG.times<0,:),1);
        icaacttmp = icaacttmp-repmat(offset,size(icaacttmp,1),1);
        mean_icacmp = double(mean(icaacttmp,2));
        latidx = dsearchn(EEG.times',lats')';
        [maxpeakvalue,maxpeak] = max(mean_icacmp(latidx(1):latidx(2)));
        [minpeakvalue,minpeak] = min(mean_icacmp(latidx(1):latidx(2)));
        maxpeak = maxpeak+latidx(1)-1;
        minpeak = minpeak+latidx(1)-1;
        %figure;plot(mean_icacmp); hold on; scatter(maxpeak,maxpeakvalue);
        h1 = kstest(icaacttmp(maxpeak,:));
        h2 = kstest(icaacttmp(minpeak,:));
        maxcon = length(find(icaacttmp(maxpeak,:)>0))/size(icaacttmp,2);
        mincon = length(find(icaacttmp(minpeak,:)<0))/size(icaacttmp,2);
        if maxcon > varargin{2} || mincon > varargin{2} 
            EEG.reject.gcompreject(1,i) = 0;
        else
            EEG.reject.gcompreject(1,i) = 1;
        end    
    end
elseif option==3 && ~isempty(varargin)
    lats = varargin{1}*1000;
    for i = 1:size(EEG.icawinv,2)
        icaacttmp  = squeeze(eeg_getdatact(EEG, 'component', i)); % provides dp x trials
        mean_icacmp = double(mean(icaacttmp,2));
        latidx = dsearchn(EEG.times',lats')';
        stdratio = std(mean_icacmp(latidx(1):latidx(2)))/std(mean_icacmp(EEG.times<0));
        if stdratio > varargin{2} && i<=varargin{3}
            EEG.reject.gcompreject(1,i) = 0;
        else
            EEG.reject.gcompreject(1,i) = 1;
        end    
    end
elseif option==4 && ~isempty(varargin)
    lats = varargin{1}*1000;
    for i = 1:size(EEG.icawinv,2)
        icaacttmp  = squeeze(eeg_getdatact(EEG, 'component', i)); % provides dp x trials
        offset     = nan_mean(icaacttmp(EEG.times<0,:),1);
        icaacttmp = icaacttmp-repmat(offset,size(icaacttmp,1),1);
        mean_icacmp = double(mean(icaacttmp,2));
        latidx = dsearchn(EEG.times',lats')';
        stdtrials = std(icaacttmp');
        nstd_base = mean_icacmp(EEG.times<0)./stdtrials(EEG.times<0)';
        nstd_peak = mean_icacmp(latidx(1):latidx(2))./stdtrials(latidx(1):latidx(2))';
        if (mean(abs(nstd_peak)) > mean(abs(nstd_base))) && i<=varargin{3}
            EEG.reject.gcompreject(1,i) = 0;
        else
            EEG.reject.gcompreject(1,i) = 1;
        end
    end
elseif option==5 && ~isempty(varargin)
    lats = varargin{1}*1000;
    for i = 1:size(EEG.icawinv,2)
        icaacttmp  = squeeze(eeg_getdatact(EEG, 'component', i)); % provides dp x trials
        mean_icacmp = double(mean(icaacttmp,2));
        latidx = dsearchn(EEG.times',lats')';
        rho_base = corr(icaacttmp(EEG.times<0,:),mean_icacmp(EEG.times<0),'type','Spearman');
        rho_peak = corr(icaacttmp(latidx(1):latidx(2),:),mean_icacmp(latidx(1):latidx(2)),'type','Spearman');
        if (mean(rho_peak) > varargin{2}*mean(rho_base)) && i<=varargin{3}
            EEG.reject.gcompreject(1,i) = 0;
        else
            EEG.reject.gcompreject(1,i) = 1;
        end
    end
elseif option==6 && ~isempty(varargin)
    lats = varargin{1}*1000;
    base = varargin{4}*1000;
    sd_out = varargin{5};
    num_out = varargin{6};
    lineartrend = varargin{7};
    peakidx = dsearchn(EEG.times',lats')';
    baseidx = dsearchn(EEG.times',base')';
    if isempty(varargin{3})
        max_comp=size(EEG.icawinv,2);
    else
        max_comp=varargin{3};
    end
    for i = 1:size(EEG.icawinv,2)
        icaacttmp  = squeeze(eeg_getdatact(EEG, 'component', i)); % provides dp x trials
        offset     = nan_mean(icaacttmp(baseidx,:),1); % offset by peak period to look at the baseline
        icaacttmp = icaacttmp-repmat(offset,size(icaacttmp,1),1);
        mean_icacmp = double(mean(icaacttmp,2));
        icawinv = EEG.icawinv(:,i);
        sdchan = std(icawinv);
        peakrange = max(mean_icacmp(peakidx(1):peakidx(2)))-min(mean_icacmp(peakidx(1):peakidx(2)));
        baserange = max(mean_icacmp(baseidx(1):baseidx(2)))-min(mean_icacmp(baseidx(1):baseidx(2)));
        peakbase_ratio = peakrange/baserange;
        
        std_icacmp = double(std(icaacttmp'));
        std_basediff = mean(std_icacmp(baseidx(1):baseidx(2)))/mean(double(std(icaacttmp(baseidx(1):baseidx(2))))); %max(std_icacmp(baseidx(1):baseidx(2)))-min(std_icacmp(baseidx(1):baseidx(2)));
        std_peakdiff = mean(std_icacmp(peakidx(1):peakidx(2)))/mean(double(std(icaacttmp(peakidx(1):peakidx(2))))); %max(std_icacmp(peakidx(1):peakidx(2)))-min(std_icacmp(peakidx(1):peakidx(2)));
        rho_base = length(find((nonzeros(tril(corr(icaacttmp(baseidx(1):baseidx(2),:),'type','Spearman'),-1)))>0.7));
        rho_peak = length(find((nonzeros(tril(corr(icaacttmp(peakidx(1):peakidx(2),:),'type','Spearman'),-1)))>0.7));
        
        compproj = EEG.icawinv(:,i)*eeg_getdatact(EEG, 'component', i, 'reshape', '2d');
        compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
        gfpcomp = squeeze(std(mean(compproj, 3),1));
        gfp_basediff = max(gfpcomp(baseidx(1):baseidx(2)))-min(gfpcomp(baseidx(1):baseidx(2)));
        gfp_peakdiff = max(gfpcomp(peakidx(1):peakidx(2)))-min(gfpcomp(peakidx(1):peakidx(2)));
        
        std_peakbase_ratio = std_peakdiff/std_basediff;
        gfp_peakbase_ratio = gfp_peakdiff/gfp_basediff;
        rho_peakbase_ratio = rho_peak/rho_base;
        
        temp(i,1) = peakbase_ratio;
        temp(i,2) = std_peakbase_ratio;
        temp(i,3) = gfp_peakbase_ratio;
        temp(i,4) = rho_peakbase_ratio;
        temp(i,5) = length(find(icawinv>sdchan*sd_out));
        
        N = length(mean_icacmp);
        a = (1:N)/N;
        r = corr(mean_icacmp,a');
        temp(i,6)= r^2>lineartrend;
        
        %for li = 1:length(lats)/2
        %    latidx = dsearchn(EEG.times',lats((li-1)*2+1:(li-1)*2+2)')';
        %    baseidx = dsearchn(EEG.times',base((li-1)*2+1:(li-1)*2+2)')';
        %    stdratio = std(mean_icacmp(latidx(1):latidx(2)))/std(mean_icacmp(baseidx(1):baseidx(2)));
        %    stdtrials = std(icaacttmp');
        %    nstd_base = mean_icacmp(baseidx(1):baseidx(2))./stdtrials(baseidx(1):baseidx(2))';
        %    nstd_peak = mean_icacmp(latidx(1):latidx(2))./stdtrials(latidx(1):latidx(2))';
        %    rho_base = corr(icaacttmp(baseidx(1):baseidx(2),:),mean_icacmp(baseidx(1):baseidx(2)),'type','Spearman');
        %    rho_peak = corr(icaacttmp(latidx(1):latidx(2),:),mean_icacmp(latidx(1):latidx(2)),'type','Spearman');
        %    if (mean(rho_peak) > varargin{2}*mean(rho_base)) && (mean(abs(nstd_peak)) > varargin{2}*mean(abs(nstd_base))) && stdratio > varargin{2} && i<=varargin{3}
        %        EEG.reject.gcompreject(1,i) = 0;
        %    %else
        %    %    EEG.reject.gcompreject(1,i) = 1;
        %    end
        %end
        %icawinv = EEG.icawinv(:,i);
        %sd = std(icawinv);
        %temp(i,1) = length(find(icawinv>sd*sd_out));
        if (temp(i,5)<=num_out) && i<=max_comp
            EEG.reject.gcompreject(1,i) = 0;
        end
    end
    detrendcomps=temp(:,6);

    
elseif option==7 && ~isempty(varargin)
    sd_out = varargin{5};
    num_out = varargin{6};
    EEG.reject.gcompreject = zeros(1,length(EEG.icachansind));
    for i = 1:length(EEG.icachansind)
        icawinv = EEG.icawinv(:,i);
        sd = std(icawinv);
        temp(i,1) = length(find(icawinv>sd*sd_out));
        if length(find(icawinv>sd*sd_out))>num_out
            EEG.reject.gcompreject(1,i) = 1;
        end
    end
end