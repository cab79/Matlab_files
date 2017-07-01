clear all
close all

grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';%Exp1 left v righ
%grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';%Exp2
cd(filepath);
run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
subjects = subjlists(grplist);
load('timefreq_limits_ERP_evoked_FNUM_CNUM');
%load cov_RT_HL-HR-PL-PR;

%select data
statmode = 'subj';
grplistselect = 1:4;
subjinfo = subjects(grplistselect);
condlist = sublist_grp(grplistselect);
condcont = [1 1 1 1]; % contrast. All '1' mean common condition values will be collapsed. '1 -1' will subtract common condition values.
%cov = cov(grplistselect);

%define latencies
latency =  [0 0.8];%timefreq_limits.limits_all{:};%
peakdef =  [1 1];% defines which peak the latencies refer to. timefreq_limits.bins; %

%set parameters
alpha = 0.05;
numrand = 1000; 
ttesttail = 0;
testgfp = 'on'; gfpbasecorrect=0;
singlesource = 'off';
testmean = 'off';
testlat = 'off';
timeshift =0;

stat = FTstats(statmode,subjinfo,condlist,condcont,latency,cov,filepath,'alpha',alpha,'numrand',numrand,'ttesttail',ttesttail,'testgfp',testgfp,...
    'singlesource',singlesource,'testmean',testmean,'testlat',testlat,'timeshift',timeshift,'peakdef',peakdef);

plotclusters(stat);

if iscell(latency)
    clusidx = stat.posclusterslabelmat>=1;
    latind = [latency{:}];
    times = -0.2:0.004:0.796;
    times(latind(clusidx))
end
