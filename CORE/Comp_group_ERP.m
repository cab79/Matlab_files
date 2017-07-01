clear all
close all

%% Group analysis - data and conditions
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);
files = dir('*4_conds_ALLEEG.mat');
files_ana = 1:length(files)-1;
trials_ana = 1; % propotion of trials to analyse

%right hand, all probs
name = 'R_allprob';
datadd = [5 6 13 14 21 22]; % change
datsub = [7 8 15 16 23 24]; % no change

%left hand, all probs
name = 'L_allprob';
datadd = [1 2 9 10 17 18];
datsub = [3 4 11 12 19 20];


%left hand, low prob
datadd = [1 2];
datsub = [3 4];
%left hand, med prob
datadd = [9 10];
datsub = [11 12];
%left hand, equal prob
datadd = [17 18];
datsub = [19 20];


%right hand, low prob
datadd = [5 6];
datsub = [7 8];
%right hand, med prob
datadd = [13 14];
datsub = [15 16];
%right hand, equal prob
datadd = [21 22];
datsub = [23 24];

%left hand, low&med probs, 1 DC
datadd = [1 9];
datsub = [3 11];
%right hand, low&med probs, 1 DC
datadd = [5 13]; % change
datsub = [7 15]; % no change
%left hand, low&med probs, 3 DC
datadd = [2 10];
datsub = [4 12];
%right hand, low&med probs, 3 DC
datadd = [6 14]; % change
datsub = [8 16]; % no change

%left hand, med probs, 3 DC
datadd = [10];
datsub = [12];
%right hand, med probs, 3 DC
datadd = [14]; % change
datsub = [16]; % no change
    
    %left hand, low prob, 1 change
    datadd = [1];
    datsub = [3];
    
    %left hand, low prob, 3 change
    datadd = [2];
    datsub = [4];
    
    
    %right hand, low prob, 1 change
    datadd = [5];
    datsub = [7];
    
    %right hand, low prob, 3 change
    datadd = [6];
    datsub = [8];


%% Analysis / plotting
GRPEEG=struct;
numtrials = [];
for f =files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    load(fullfile(filepath,files(f).name));
    
    totadd = 0;
    for ad = 1:length(datadd)
        EEG = ALLEEG(datadd(ad));
        no_trials_add = ceil(length(EEG.epoch)*trials_ana);
        totadd = totadd+no_trials_add;
        selecttrials = 1:no_trials_add;
        ALLEEG(datadd(ad)) = pop_select(EEG,'trial',selecttrials);
    end
    
    totsub = 0;
    for ad = 1:length(datsub)
        EEG = ALLEEG(datsub(ad));
        no_trials_sub = ceil(length(EEG.epoch)*trials_ana);
        totsub = totsub+no_trials_sub;
        selecttrials = 1:no_trials_sub;
        ALLEEG(datsub(ad)) = pop_select(EEG,'trial',selecttrials);
    end
    
    numtrials(f,:)=[totadd totsub];
    
    if length(datadd)>1; 
        addEEG = pop_mergeset(ALLEEG, datadd, 1);
    else
        addEEG = ALLEEG(datadd);
    end;
    if length(datsub)>1; 
        subEEG = pop_mergeset(ALLEEG, datsub, 1); 
    else
        subEEG = ALLEEG(datsub);
    end;
    
    if f==files_ana(1)
        GRPEEG = addEEG;
    else 
        GRPEEG(f) = addEEG;
    end
    GRPEEG(f+length(files_ana)) = subEEG;
end

[erp1 erp2 erpsub, tms, pvalues, tvalues, erp1_sd, erp2_sd] = pop_comperpCAB(GRPEEG, 1, 1:length(files_ana), length(files_ana)+1:length(files_ana)*2,'addavg','on','subavg','on','diffavg','on','alpha',0.05);

erp_plot = {erp1 erp2};
erp_plotvar = {erp1_sd erp2_sd};
plotnames = {'no change' 'change'};
col = {'b','r'};
% plot data
ptimes = find(((EEG.times >= -100).*(EEG.times <= 200))==1);
ele=[34 70];
sigtimes=[];%245:274;%find(stat.posclusterslabelmat); % from ft_ClusterStatDesign
for e=1:length(ele)
    figure
    hold all
    st=EEG.times(sigtimes);
    ts=EEG.times(ptimes);
    fill([st, fliplr(st)], [ones(1,length(st))*2, fliplr(ones(1,length(st))*-1.5)], ...
    [0.5 0.5 0.5], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    for ii = 1:length(erp_plot); 
        peakdata = erp_plot{ii};
        vardata = erp_plotvar{ii};
        ERP = mean(peakdata(ele(e),:),1);
        VAR = mean(vardata(ele(e),:),1);
        nsub = length(files_ana);
        SEM = VAR/sqrt(nsub);               % Standard Error
        tscore = tinv(0.025,nsub-1);      % T-Score
        CI = tscore*SEM;                      % Confidence Intervals
        upper = ERP(ptimes)+CI(ptimes);
        lower = ERP(ptimes)-CI(ptimes);
        fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
        col{ii}, 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
        plot(ts,ERP(ptimes),col{ii}); 
    end
    %legend('Healthy controls', 'healthy unaff', 'CRPS patients','patient unaff')
    %legend('Healthy controls', 'CRPS patients','patient unaff');
    ylabel('Amplitude, uV') % label for y axis
    xlabel('Time (ms)') % label for x axis
    set(gca,'FontSize',15);
    set(findall(gcf,'type','text'),'FontSize',15);
end

%sname = ['grpdata_' name];
%save(sname,'GRPEEG','-v7.3');


% plot grand average
plotnames = {'no change' 'change' 'subtracted'};
all_plot = {erp1 erp2 erpsub};
for ii = 1:length(all_plot);
    plotdata = all_plot{ii};
    if ~strcmp(plotnames{ii},'subtracted')
        [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
        [~, maxmaxidx] = max(maxval);
        plottime = EEG.times(maxidx(maxmaxidx));
        if plottime == EEG.times(end)
            plottime = EEG.times(end-1);
        end
    end
    figure('Name',plotnames{ii});
    timtopo(plotdata,EEG.chanlocs,...
        'limits',[EEG.times(1) EEG.times(end)],...
        'plottimes',plottime,'maplimits',[-1 1]);
    set(gcf,'Color','white');
end