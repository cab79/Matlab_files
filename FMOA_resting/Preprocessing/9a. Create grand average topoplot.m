load gamma1;
load gamma2;
load gamma3;
load beta;
load alpha1;
load alpha2;
load theta;
load delta;

gamma1 = mean(gamma1,2);
gamma1 = mean(gamma1,3);
gamma2 = mean(gamma2,2);
gamma2 = mean(gamma2,3);
gamma3 = mean(gamma3,2);
gamma3 = mean(gamma3,3);
beta = mean(beta,2);
beta = mean(beta,3);
alpha1 = mean(alpha1,2);
alpha1 = mean(alpha1,3);
alpha2 = mean(alpha2,2);
alpha2 = mean(alpha2,3);
theta = mean(theta,2);
theta = mean(theta,3);
delta = mean(delta,2);
delta = mean(delta,3);

topoplot(gamma1([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('gamma1'); %caxis([-0.1 0.1])

topoplot(gamma2([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('gamma2'); %caxis([-0.1 0.1])

topoplot(gamma3([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('gamma3'); %caxis([-0.1 0.1])

topoplot(beta([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('beta'); %caxis([-0.1 0.1])

topoplot(alpha1([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('alpha1'); %caxis([-0.1 0.1])

topoplot(alpha2([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('alpha2'); %caxis([-0.1 0.1])

topoplot(theta([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('theta'); %caxis([-0.1 0.1])

topoplot(delta([1:2 4:30 33:64], 1), 'C:\Documents and Settings\mdrabek\Desktop\EEG analysis Matlab\chan61.locs'), colorbar, title('delta'); %caxis([-0.1 0.1])
