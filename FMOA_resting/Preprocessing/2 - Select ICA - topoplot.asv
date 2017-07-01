
% get plot of head with source of component
clear all
subject = 'OA17';
load ([subject '_data_epochs.mat']);
load(['b_' subject]);
act = b * total_data; 

for i = 1:40
act = b * total_data;
ib = pinv(b);
iact= ib(:, i) * act(i,:);  
topoplot(mean(iact,2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs'); caxis([-0.1 0.1]), colorbar, title(num2str(i))
pause
end