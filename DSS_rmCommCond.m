function [newdata] = DSS_rmCommCond(data,iconds,ord,ncomp,plottoggle)
% removes common evoked activity from a pair or multiple conditions. Ideally use after DSS has already been applied to de-noise ERP for each condition. 

% if data is a cell array of eeg data matrices from different conditions,
% iconds should be empty and will created in this function

% iconds: array of condition numbers, length of trials. Zeros will prevent
% those trials from being processed.

if iscell(data)
    if ~isempty(iconds)
        error('cannot differentiate conditions within each of multiple data matrices');
    end
    ntrials=cell(1,length(data));
    datacat=[];
    for d = 1:length(data)
        data{d} = permute(data{d},ord);
        if~isempty(data{d}); [nsamps,nchans,ntrials{d}] = size(data{d});end;
        datacat = cat(3,datacat,data{d});
        iconds = [iconds d*ones(1,ntrials{d})];
    end
else
    datacat = permute(data,ord);
    [nsamps,nchans,ntrials] = size(datacat);
end

newdata=data;
newdatacat=datacat;
%if isempty(iconds); iconds = ones(1,ntrials);end;
numcond = sort(unique(iconds));

% apply DSS to clean them
c1=zeros(nchans);
[c1,nc1] = nt_cov(mean(datacat(:,:,iconds>0),3));
c1=c1/nc1;
p1=cell(1);
z=cell(1);
for iCondition=numcond
    c0=zeros(nchans); 
    [c0,nc0] = nt_cov(datacat(:,:,iconds==iCondition));
    c0=c0/nc0;
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    p1{iCondition}=pwr1./pwr0; % score, proportional to power ratio of ERP to full data
    z{iCondition}=nt_mmat(datacat(:,:,iconds==iCondition),todss);
    %newdatacat(:,:,iconds==iCondition) = nt_mmat(z(:,:,iconds==iCondition),todss');
    newdatacat(:,:,iconds==iCondition)=nt_tsr(datacat(:,:,iconds==iCondition),z{iCondition}); % regress them out
end

% plot results
if strcmp(plottoggle,'on')
    for iCondition=numcond
        figure
        % plot bias score
        clf; set(gcf,'color', [1 1 1]);
        subplot 221; 
        plot(p1{iCondition}, '.-'); xlabel('component'); ylabel('score'); title('DSS score');
        subplot 222; 
        plot(mean(datacat(:,:,iconds==iCondition),3)); title(['data cond' num2str(iCondition)]);
        subplot 223;
        nt_bsplot(mean(z{iCondition}(:,1:ncomp,:),2)); title('removed'); 
        subplot 224;
        plot(mean(newdatacat(:,:,iconds==iCondition),3)); title('reconstructed'); 
    end
end

if iscell(data)
    tend=0;
    for d = 1:length(data)
        tstart = tend+1;
        tend = tend+ntrials{d};
        newdata{d} = newdatacat(:,:,tstart:end);
        newdata{d} = permute(newdata{d},ord);
    end
else
    newdata = permute(datacat,ord);
end
