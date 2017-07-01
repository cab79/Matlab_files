clear all

fnames = {'delta','theta','alpha','beta','gamma1','gamma2'};
fname2 = '_ANOVAele_stats.mat';
chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';

for f = 1:length(fnames)
    fname = [fnames{f} fname2];
    load(fname);
    
    subplot(1,3,1); 
    [handle,Zi,grid] = topoplot(FvPv(:,1,1), chan_locs, 'maplimits','absmax'); 
    colorbar
    title('time');
    
    subplot(1,3,2); 
    [handle,Zi,grid] = topoplot(FvPv(:,1,2), chan_locs, 'maplimits','absmax'); 
    colorbar
    title('scan');
    
    subplot(1,3,3); 
    [handle,Zi,grid] = topoplot(FvPv(:,1,3), chan_locs, 'maplimits','absmax'); 
    colorbar
    title('time*scan');
    
    mtit(fnames{f});
    pause
end