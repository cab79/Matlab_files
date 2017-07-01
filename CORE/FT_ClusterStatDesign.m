clear all
close all

%% Data info
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);
load('C:\Data\CORE\Preprocessed_100Hz\groupsublist.mat');
subjlists={};
sub_prefix = {'CORE'}; 
file_suffix = '_4_conds_ALLEEG.mat'; 
grplist = [1]; 
%sublist_side = {'LR'}; 
%expts = {'4'}; 
%hand_nme = {'L','R'};

%% Conditions to compare
%left hand, low prob
%cond1 = [1 2]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%cond2 = [3 4]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%condlist = {'Ch' 'NoCh'};
%chan_select = 'E93';

%right hand, low prob
cond1 = [5 6]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
cond2 = [7 8]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
condlist = {'Ch' 'NoCh'};
chan_select = 'E42';

%left hand, equal prob
%cond1 = [17 18]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%cond2 = [19 20]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%condlist = {'Ch' 'NoCh'};
%chan_select = 'E93';

%right hand, equal prob
%cond1 = [21 22]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%cond2 = [23 24]; % ROWS : difference waveform (ERP) tested; COLUMNS: averaged
%condlist = {'Ch' 'NoCh'};
%chan_select = 'E42';


%% Stats parameters
trials_ana = 1; % propotion of trials to analyse
numERPcomp=[];
covariate = [];
%---select tests and data---%
%statmode = 'subj_corr'; grplistselect = [1;2;3;4]; condcont = [1]; numtests = 4; 
%statmode = 'corr'; condlist = sublist_grp(i); covariate = cov(:,i); 
%statmode = 'subj'; grplistselect = [1 2 3 4]; condlist = sublist_grp(grplistselect); numtests = 1; condcont = [1 -1 1 -1]; 
statmode = 'cond'; numtests = 1; condcont = [1 -1]; use_etype =[];% grplistselect = one row per numtest. condcont = ones unless subtraction required
%statmode = 'condtrial'; grplistselect = [1]; numtests = 1; condlist = {'2' '10'}; condcont = [1 -1];

%% Create group / subject list structure
for g = 1:size(group,2)
    subgrp={};
    s2=0;
    for s = 1:size(group,1)
        if isnan(group(s,g))
            continue
        end
        s2 = s2+1;
        subj = num2str(group(s,g));
        if length(subj)<2
            subj = ['00' subj];
        elseif length(subj)<3
            subj = ['0' subj];
        end
        subgrp{s2,1} = [sub_prefix{g} subj];
    end
    subjlists{g,1} = subgrp;
end
subjects = subjlists(grplist);



%% Create data structure
subjinfo=cell(1,1);
numtrials = [];

if size(cond1,1)==2
    condcont = [1 -1];
end

for g = 1:length(subjects)
    for s = 1:length(subjects{g,1}) 
        subj = subjects{g,1}{s,1};
        load(fullfile(filepath,[subj file_suffix]));
    
        %--- to test a subset of trials only ---%
        %totadd = 0;
        %for ad = 1:length(cond1)
        %    EEG = ALLEEG(cond1(ad));
        %    no_trials_add = ceil(length(EEG.epoch)*trials_ana);
        %    totadd = totadd+no_trials_add;
        %    selecttrials = 1:no_trials_add;
        %    ALLEEG(cond1(ad)) = pop_select(EEG,'trial',selecttrials);
        %end

        %totsub = 0;
        %for ad = 1:length(cond2)
        %    EEG = ALLEEG(cond2(ad));
        %    no_trials_sub = ceil(length(EEG.epoch)*trials_ana);
        %    totsub = totsub+no_trials_sub;
        %    selecttrials = 1:no_trials_sub;
        %    ALLEEG(cond2(ad)) = pop_select(EEG,'trial',selecttrials);
        %end
        %numtrials(g,s,:)=[totadd totsub];

        for diff = 1:max(size(cond1,1),size(cond2,1))
            if size(cond1,2)>1; 
                cond1EEG = pop_mergeset(ALLEEG, cond1(diff,:), 1);
            else
                cond1EEG = ALLEEG(cond1(diff,:));
            end;
            if size(cond2,2)>1; 
                cond2EEG = pop_mergeset(ALLEEG, cond2(diff,:), 1); 
            else
                cond2EEG = ALLEEG(cond2(diff,:));
            end;

            subjinfo{g,1}{s,1}{1,1}{diff,1} = cond1EEG;
            subjinfo{g,1}{s,1}{2,1}{diff,1} = cond2EEG;
        end
    end
end

%% Run FT stats
statall = cell(1,numtests);
for i = 1:numtests
    
    if strcmp(statmode,'corr') || strcmp(statmode,'subj_corr') 
        condlist = sublist_grp(i); condlist = [condlist 'cov']; 
        covariate = cov(:,i); 
    %elseif strcmp(statmode,'cond')
    %    condlist = sublist_side(grplistselect(i,:));
    end

    %define latencies
    latency =  [0 0.1];%timefreq_limits.limits_all;%
    peakdef =  {[1 1]};%timefreq_limits.bins;%            % defines which peak the latencies refer to. MUST BE CELL ARRAY
    frequency = [];%4:2:20;%[10]; % empty if ERP, performs freq analysis if specified

    %set parameters
    alpha = 0.05;
    numrand = 1000; 
    ttesttail = 0;
    test_gfporchan = chan_select; gfpbasecorrect=0; % option: gfp, chan number, off
    singlesource = 'off';
    testmean = 'off';
    testlat = 'off';
    timeshift =0;

    stat = FTstats(statmode,subjinfo,condlist,condcont,latency,frequency,covariate,filepath,'alpha',alpha,'numrand',numrand,'ttesttail',ttesttail,'test_gfporchan',test_gfporchan,...
        'singlesource',singlesource,'testmean',testmean,'testlat',testlat,'timeshift',timeshift,'peakdef',peakdef,'numERPcomp',numERPcomp);
    

    statall{i} = stat;

    %if iscell(latency)
    %    clusidx = stat.posclusterslabelmat>=1;
    %    latind = [latency{:}];
    %    times = -0.2:0.004:0.796;
    %    times(latind(clusidx))
    %end

end

%% Plot results
close all
for i = numtests
    stat = statall{1,i};
    if ~strcmp(test_gfporchan,'off');
        stat.cfg.alpha = 0.05;
        plotclusters(stat);
    else
        cfg = [];
        cfg.alpha  = 0.05;
        cfg.zlim = [-6 6]; %Tvalues
        cfg.elecfile = 'FT_layout.mat';
        %cfg.parameter = 'diffcond';
        
        if ~isempty(frequency)
            for f = 1:size(stat.diffcond.avg,2)
                statf = stat;
                statf.posclusterslabelmat = squeeze(statf.posclusterslabelmat(:,f,:));
                statf.negclusterslabelmat = squeeze(statf.negclusterslabelmat(:,f,:));
                statf.prob = squeeze(statf.prob(:,f,:));
                statf.cirange = squeeze(statf.cirange(:,f,:));
                statf.mask = squeeze(statf.mask(:,f,:));
                statf.stat = squeeze(statf.stat(:,f,:));
                statf.ref = squeeze(statf.ref(:,f,:));
                %statf.diffcond = squeeze(statf.diffcond.avg(:,f,:));
                statf = rmfield(statf,'freq');
                ft_clusterplot(cfg,statf)
            end
        else
            %stat.diffcond = stat.diffcond.avg;
            ft_clusterplot(cfg,stat)
        end
    end
    if strcmp(statmode,'corr') || strcmp(statmode,'subj_corr') 
        latidx = dsearchn(stat.diffcond.time',latency')';
        data = stat.diffcond.individual(:,:,latidx(1):latidx(2));
        numcls = length(find([stat.posclusters.prob]<alpha));
        for ns = 1:numcls
            dep = mean(data(:,stat.posclusterslabelmat==ns),2);
            figure
            scatter(cov(:,i),dep);
            [r p] = corr(cov(:,i),dep,'type', 'Spearman')
        end
    end
end