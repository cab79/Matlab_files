clear all
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);
files = dir('*_conds_ALLEEG.mat');
trials_ana = 1; % propotion of trials to analyse

%% plot conditions
for f =57%:length(files)
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    load(fullfile(filepath,files(f).name));
    %left hand, all probs
    datadd = [1 2 9 10 17 18];
    datsub = [3 4 11 12 19 20];
    
    %left hand, low prob
    datadd = [1 2];
    datsub = [3 4];
    
    %left hand, low prob, 1 change
    datadd = [1];
    datsub = [3];
    
    %left hand, low prob, 3 change
    datadd = [2];
    datsub = [4];
    
    %left hand, equal prob
    datadd = [17 18];
    datsub = [19 20];
    
    %right hand, all probs
    datadd = [5 6 13 14 21 22];
    datsub = [7 8 15 16 23 24];
    
    %right hand, low prob
    datadd = [5 6];
    datsub = [7 8];
    
    %right hand, equal prob
    datadd = [21 22];
    datsub = [23 24];
    
    %right hand, low prob, 1 change
    datadd = [5];
    datsub = [7];
    
    %right hand, low prob, 3 change
    datadd = [6];
    datsub = [8];
    
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
    
    numtrials = [];
    numtrials(f,:)=[totadd totsub];
    
    [erp1 erp2 erpsub] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on','ylim',[-3 3]);
end

%% plot digits
%left hand, 3 digit change
datadd = [1];
datsub = [4];

%right hand, 3 digit change
datadd = [5];
datsub = [8];

%left hand, 1 digit change
datadd = [2];
datsub = [3];

%right hand, 2 digit change
datadd = [6];
datsub = [7];

[erp1 erp2 erpsub] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on');

close all
for digit = 1:4
    figure
    pop_timtopo(ALLEEG(digit), [-100 298], [60 90 110], ['digit_' num2str(digit)], 'maplimits', [-2 2]);
end

%% Group analysis
datadd = [1];
datsub = [2];

[erp1 erp2 erpsub sig] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on','alpha',0.05);
