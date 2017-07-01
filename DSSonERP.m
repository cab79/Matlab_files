function newdata = DSSonERP(data,ord,iconds,keep1,keep2,plottoggle)
% iconds: array of condition numbers, length of trials

% Same as example1, but the data now include multiple conditions.
% We look for the linear combination that maximizes repeatability jointly
% for all conditions.  Data are in a cell array of matrices of dimensions 
% time*channels*trials
%
% Uses nt_dss0().


data = permute(data,ord);
[nsamps,nchans,ntrials] = size(data);
if isempty(iconds); iconds = ones(1,ntrials);end;
numcond = sort(unique(iconds));
newdata=data;

% apply DSS to clean them
c0=zeros(nchans); c1=zeros(nchans);
c0_all =c0;
c1_all =c1;
for iCondition=numcond
    [c0 nc0] = nt_cov(data(:,:,iconds==iCondition));
    c0=c0/nc0;
    c0_all = c0_all+c0;
    [c1 nc1] = nt_cov(mean(data(:,:,iconds==iCondition),3));
    c1=c1/nc1;
    c1_all = c1_all+c1;
end
[todss,pwr0,pwr1]=nt_dss0(c0_all,c1_all,keep1,keep2);
p1=pwr1./pwr0; % score, proportional to power ratio of ERP to full data
for iCondition=numcond
    z(:,:,iconds==iCondition)=nt_mmat(data(:,:,iconds==iCondition),todss);
    newdata(:,:,iconds==iCondition) = nt_mmat(z(:,:,iconds==iCondition),todss');
end

% plot results
if strcmp(plottoggle,'on')
    % plot bias score
    figure; clf; set(gcf,'color', [1 1 1]);
    plot(p1, '.-'); xlabel('component'); ylabel('score'); title('DSS score');

    for iCondition=numcond
        figure
        subplot 131; 
        plot(mean(data(:,:,iconds==iCondition),3)); title(['data cond' num2str(iCondition)]);
        subplot 132;
        nt_bsplot(mean(z(:,1:keep1,iconds==iCondition),2)); title('recovered'); 
        subplot 133;
        plot(mean(newdata(:,:,iconds==iCondition),3)); title('reconstructed'); 
    end
end

newdata = newdata(:,:,iconds>0); % only keep trials specified by iconds
newdata = ipermute(newdata,ord);
