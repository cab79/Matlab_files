%% plot conditions for each subject - TO IDENITFY THOSE SUBJECTS WHO REQUIRE FURTHER CLEANING / REJECTING
% Requires ALLEEG files saved from preprocess part 2

clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*_ALLEEG.mat');

% To use this script, just type one number at a time here (index of "files")
files_ana = [1]; % values from 1 to length(files)

for f = files_ana
    
    % load by index in list
    load(fullfile(filepath,files(files_ana(f)).name))
   
    % choose a maximum of two condition to plot on one chart. Best to
    % choose those conditions most likely to be noisy. Will need to run a
    % second time to compare a different two conditions. This example has 8
    % conditions and the first two will be plotted for now.
    datadd = [1];
    datsub = [2];
    
    %datadd = [3];
    %datsub = [4];
    
    %datadd = [5];
    %datsub = [6];
    
    %datadd = [7];
    %datsub = [8];
    
    %close all
    [erp1 erp2 erpsub] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on');
    %pause
end


%for cond = 1:2
%    figure
%    pop_timtopo(ALLEEG(cond), [-3000 1000], [60 90 110], ['cond_' num2str(cond)], 'maplimits', [-2 2]);
%end
