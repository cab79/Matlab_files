function S=MultiOutliers(S,data);
% Between-subject analysis: assumes data is a cell array of Fieldtrip data structures containing .avg
% field
% Within-subject analysis: assumes data is a 3D double array of
% chan_time_trials

%restoredefaultpath
addpath('C:\Data\Matlab\Multivariate_Outliers')

%X = %N*d. N - number of samples (subjects (e.g. 120) * trials (e.g. 100) = 12,000), d is dimensionality (channels (e.g. 59) * data points per trial (e.g. 6000/4 = 1500) = 100,000)
% not possible, so can have no more than a selection of data points:
% Nsub * Ntrial / Nchan.
% can either run multiple times over time windows, or just use random
% selection of time points.
%outliers_demo_v2(12000,64,0.05)

% for one data point at a time:
% Subs (e.g. 30) * conds (e.g. 2) / chans (e.g. 59) need to be greater than
% 1.

% Over trials within-subject:
% Trials (e.g. 400) / chans (e.g. 92)

% create data array X
if iscell(data) && isfield(data{1},'avg') % SUBJECTS
    dattype='subjects'
    for f = 1:length(data)
        dat=data{f}.avg(S.(S.func).inclchan,:);
        for e = 1:size(dat,1)
            dat(e,:) = smooth(squeeze(dat(e,:)),30,'lowess');
        end
        %dat2 = ft_preproc_bandpassfilter(dat, 250, [0 20]);
        X{f}=dat;
    end
    X = cat(3,X{:});
elseif isnumeric(data) && ndims(data)==3 % TRIALS
    dattype='trials'
    X=data;
end

% cycle through for each time point
for i = 1:size(X,2)
    Xi = squeeze(X(:,i,:))';
    %close all
    [mu,var,RD,chi_crt]=DetectMultVarOutliers(Xi,[],[],0);
    S.(S.func).multout.mu(i,:) = mu;
    S.(S.func).multout.var(i,:,:) = var;
    S.(S.func).multout.RD(i,:) = RD';
    S.(S.func).multout.chi_crt(i,:) = chi_crt;
    disp(['sample ' num2str(i) ' of ' num2str(size(X,2))])
end

% plots
figure; plot(S.(S.func).multout.mu);
title('robust mean')

hf=figure('color','w');
drawnow
pause(0.1)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
j_hf=get(hf,'JavaFrame'); 
j_hf.setMaximized(true);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
drawnow
pause(0.1)

RD=log10(median(S.(S.func).multout.RD,1));
chi_crt=log10(median(S.(S.func).multout.chi_crt,1));
N=length(RD);
spl_id=1:N;
h_a=plot(spl_id,RD','ok','MarkerSize',6,'MarkerFaceColor','w'); % basis subset
hold on
%h_b=plot(spl_id(~id_in),RD(~id_in)','ok','MarkerSize',6,'MarkerFaceColor','k'); % complement of the basis subset

if N<=10
    dx=0.25;
elseif N<=100
    dx=0.5;
elseif N<=250
    dx=1;
elseif N<=500
    dx=2;
elseif N<=1E3
    dx=4;
else
    dx=10;
end

h1=plot([1-dx N+dx],chi_crt(1)*[1 1],':','Color',[0.75 0 0.75],'LineWidth',1);
h2=plot([1-dx N+dx],chi_crt(2)*[1 1],':','Color',[0.75 0 0],'LineWidth',1);
h3=plot([1-dx N+dx],chi_crt(3)*[1 1],':','Color',[0 0.75 0],'LineWidth',1);
h4=plot([1-dx N+dx],chi_crt(4)*[1 1],':','Color',[0 0 0.75],'LineWidth',1);

xl=sprintf('sample id (1 to %u)',N);
xlabel(xl,'FontSize',25,'FontWeight','bold','Color','k')
ylabel('log Robust Mahalanobis Distance','FontSize',25,'FontWeight','bold','Color','k')
set(gca,'FontSize',20,'XLim',[0 N+1],'XColor','k','YColor','k')

xt=get(gca,'XTick');
xt(1)=1;
xt=unique(xt);
set(gca,'XLim',[1-dx N+dx],'XTick',xt,'TickDir','out')

h=legend([h_a h1 h2 h3 h4],...
         {'data' '$$\alpha=0.2$$' '$$\alpha=0.1$$' '$$\alpha=0.05$$' '$$\alpha=0.01$$'},...
          'Location','EastOutside',...
          'Orientation','vertical',...
          'Interpreter','latex');

p=get(h,'Position');
x=1-1.02*p(3); p(1)=x; set(h,'Position',p)

drawnow
pause(0.1)

p=get(gca,'Position'); w=0.98*(x-p(1)); p(3)=w; set(gca,'Position',p)

if strcmp(dattype,'subjects')
    % summarise per subject
    files = unique(S.(S.func).fileidx);
    subs = S.(S.func).designmat(2:end,find(strcmp(S.(S.func).designmat(1,:),'subjects')));
    [unisubs,~,subsidx] = unique(subs,'stable');
    for i = 1:length(unisubs)
        RDsub(i) = mean(RD(ismember(S.(S.func).fileidx,find(subsidx==i))));
    end
    [sortRD,ord] = sort(RDsub,'descend');
    S.(S.func).multoutlist = unisubs(ord);
    S.(S.func).multoutlist(:,2) = num2cell(sortRD');
elseif strcmp(dattype,'trials')
    [sortRD,ord] = sort(RD,'descend');
    S.(S.func).multoutlist = ord;
    S.(S.func).multoutlist(:,2) = sortRD;
end