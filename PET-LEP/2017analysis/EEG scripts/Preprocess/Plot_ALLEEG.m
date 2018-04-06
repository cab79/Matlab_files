clear all
filepath = 'C:\Data\PET-LEP\Preprocessed';
cd(filepath);
files = dir('*_ALLEEG.mat');

%% plot conditions
for f = 3%1:length(files)
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    load(fullfile(filepath,files(f).name));
   
    datadd = [1];
    datsub = [2];
    
    [erp1 erp2 erpsub] = pop_comperp(ALLEEG, 1, datadd, datsub,'addavg','on','subavg','on','diffavg','on');
    %pause
end


%close all
for cond = 1:2
    figure
    pop_timtopo(ALLEEG(cond), [-3000 1000], [100:50:600], ['cond_' num2str(cond)], 'maplimits', [-2 2]);
end
