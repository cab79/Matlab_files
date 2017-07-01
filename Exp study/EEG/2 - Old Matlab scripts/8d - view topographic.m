



grand_all=zeros(64,625);

fnames={'1_1_grand_avg.mat';
    '1_2_grand_avg.mat';
    '2_1_grand_avg.mat';
    '2_2_grand_avg.mat';
    '3_1_grand_avg.mat';
    '3_2_grand_avg.mat'};

for x=1:length(fnames);
    
        
    for n = 1
	fname=char(fnames(x));

        load(fname);
        grand_all=grand_all + grand_avg;
        
        
    end
    grand_all = grand_all/6;
fname=char(fnames(x));    
eval(['save ' 'grand_avg.mat'  ' grand_avg'])

end









topoplot(mean(grand_all([1:2 4:30 33:64], 126:188),2), 'E:\Debbie_Backup\Matlab scripts\eeglab4.08\chan.locs2'); caxis([-0.5 0.5]), colorbar

topoplot(mean(grand_all([1:2 4:30 33:64], 251:313),2), 'E:\Debbie_Backup\Matlab scripts\eeglab4.08\chan.locs2'); caxis([-3 3]), colorbar

topoplot(mean(grand_all([1:2 4:30 33:64], 375:438),2), 'E:\Debbie_Backup\Matlab scripts\eeglab4.08\chan.locs2'); caxis([-3 3]), colorbar

topoplot(mean(grand_all([1:2 4:30 33:64], 480:500),2), 'E:\Debbie_Backup\Matlab scripts\eeglab4.08\chan.locs2'); caxis([-1 1]), colorbar



