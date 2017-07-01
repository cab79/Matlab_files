clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*_ALLEEG.mat');


files_ana = [1];%1:length(files);

%% plot conditions for each subject - IDENITFY THOSE WHO REQUIRE FURTHER CLEANING
for f = 1%[10 18 24 25 34 35 39]%:length(files)
    
    %by index in list
    load(fullfile(filepath,files(files_ana(f)).name));
    
   
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
