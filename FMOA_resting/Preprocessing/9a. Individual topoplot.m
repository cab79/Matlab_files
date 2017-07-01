clear all



load gamma1_closed;
load gamma2_closed;
load gamma3_closed;
load beta_closed;
load alpha1_closed;
load alpha2_closed;
load theta_closed;
load delta_closed;

data = alpha1_closed;

%ref_ele = [1:2 4:30 33:64];
ref_ele = [1 2 33 36 44 45 42 34 38];

refdata = mean(data(ref_ele,:),1);
data = data - repmat(refdata,size(data,1),1);

for i = 1:size(data,2)
topoplot(data([1:2 4:30 33:64], i), 'H:\Data analysis files\EEG analysis Matlab\chan.locs2'), colorbar, title('alpha1'); %caxis([-0.1 0.1])
pause;
end