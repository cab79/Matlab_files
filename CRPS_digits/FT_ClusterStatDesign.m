clear all
close all

%grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';%Exp1 left v right
grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';%Exp2
cd(filepath);
run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
subjects = subjlists(grplist);
load('timefreq_limits_ERP_evoked_FNUM_CNUM');
load cov_RT_HL-HR-PL-PR;
eventtypes = {'FNUM','CNUM','ACCU'; [1:5], [0 1], [0 1]}; use_etype = [1 2];

numERPcomp=[];
covariate = [];
%---select tests and data---%
%statmode = 'subj_corr'; grplistselect = [1;2;3;4]; condcont = [1]; numtests = 4; 
%statmode = 'corr'; condlist = sublist_grp(i); covariate = cov(:,i); 
statmode = 'subj'; grplistselect = [1 2 3 4]; condlist = sublist_grp(grplistselect); numtests = 1; condcont = [1 -1 1 -1]; 
%statmode = 'cond'; grplistselect = [1 2;3 4]; numtests = 2; condcont = [1 1];
%statmode = 'condtrial'; grplistselect = [1]; numtests = 1; condlist = {'2' '10'}; condcont = [1 -1];


statall = cell(1,numtests);
for i = 1:numtests
    
    if strcmp(statmode,'corr') || strcmp(statmode,'subj_corr') 
        condlist = sublist_grp(i); condlist = [condlist 'cov']; 
        covariate = cov(:,i); 
    elseif strcmp(statmode,'cond')
        condlist = sublist_side(grplistselect(i,:));
    end
    
    subjinfo = subjects(grplistselect(i,:));

    %define latencies
    latency =  [0 0.8];%timefreq_limits.limits_all;%
    peakdef =  {[1 1]};%timefreq_limits.bins;%            % defines which peak the latencies refer to. MUST BE CELL ARRAY
    frequency = [];%4:2:20;%[10]; % empty if ERP, performs freq analysis if specified

    %set parameters
    alpha = 0.05;
    numrand = 1000; 
    ttesttail = 0;
    testgfp = 'off'; gfpbasecorrect=0;
    singlesource = 'off';
    testmean = 'off';
    testlat = 'off';
    timeshift =0;

    stat = FTstats(statmode,subjinfo,condlist,condcont,latency,frequency,covariate,filepath,'alpha',alpha,'numrand',numrand,'ttesttail',ttesttail,'testgfp',testgfp,...
        'singlesource',singlesource,'testmean',testmean,'testlat',testlat,'timeshift',timeshift,'peakdef',peakdef,'numERPcomp',numERPcomp,'eventtypes',eventtypes,'use_etype',use_etype);
    
    statall{i} = stat;


    %if iscell(latency)
    %    clusidx = stat.posclusterslabelmat>=1;
    %    latind = [latency{:}];
    %    times = -0.2:0.004:0.796;
    %    times(latind(clusidx))
    %end

end
close all
for i = numtests
    stat = statall{1,i};
    if strcmp(testgfp,'on');
        plotclusters(stat);
    else
        cfg = [];
        cfg.alpha  = 0.05;
        cfg.zlim = [-6 6]; %Tvalues
        cfg.alpha = 0.05;
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