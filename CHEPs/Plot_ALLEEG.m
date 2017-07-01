clear all
filepath = 'C:\Data\CHEPs\Preprocessed';
cd(filepath);
files = dir('*ALLEEG.mat');

%% plot conditions
for f = 1%:length(files)
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    load(fullfile(filepath,files(f).name));
    %left hand, all probs
    timtopo(mean(ALLEEG.data,3),ALLEEG.chanlocs);
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

%% Compare 
datadd = [1];
datsub = [2];

[erp1 erp2 erpsub] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on');
