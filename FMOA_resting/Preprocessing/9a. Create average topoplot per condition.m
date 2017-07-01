clear all

load gamma2;

for sub = 1:29;

gamma_open = mean(gamma2(:,sub,[1 3]),3);
gamma_closed = mean(gamma2(:,sub,[2 4]),3);

topoplot(gamma_open([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title([num2str(sub) ' open']); caxis([-10 10])
pause
topoplot(gamma_closed([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title([num2str(sub) ' closed']); caxis([-10 10])
pause
end


%% Grand avg per cond
clear all

load gamma2;

group = [1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2];

gp1 = find(group == 1);
open_gp1_pre = mean(gamma2(:,gp1,[1]),2);
open_gp1_post= mean(gamma2(:,gp1,[3]),2);
closed_gp1_pre = mean(gamma2(:,gp1,[2]),2);
closed_gp1_post = mean(gamma2(:,gp1,[4]),2);

gp2 = find(group == 2);
open_gp2_pre = mean(gamma2(:,gp2,[1]),2);
open_gp2_post = mean(gamma2(:,gp2,[3]),2);
closed_gp2_pre = mean(gamma2(:,gp2,[2]),2);
closed_gp2_post = mean(gamma2(:,gp2,[4]),2);

topoplot(open_gp1_pre([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title(['open_gp1_pre']); caxis([-20 0])
pause
topoplot(open_gp1_post([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title(['open_gp1_post']); caxis([-20 0])
pause
topoplot(open_gp2_pre([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title(['open_gp2_pre']); caxis([-20 0])
pause
topoplot(open_gp2_post([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title(['open_gp2_post']); caxis([-20 0])
pause