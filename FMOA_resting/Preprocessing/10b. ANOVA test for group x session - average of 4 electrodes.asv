clear all
% INPUT Args:
subjects = {'M1';'M2';'M4';'M6';'M7';'M9';'M10';'M12';'M14';'M15';'M16';'M17';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M28';'M29';'M30';'M32';'M33';'M35';'M36';'M37';'M38';'M40'};
NSub = length(subjects);
Ngp1 = 15;
Ngp2 = 14;
%ele = [1:2 4:30 33:64];
ele = [62 30 57 24];
results = zeros(3,3,1);

%%%%%%%%%%

subject = [1:NSub 1:NSub]; %- subject ID
time = [ones(1,NSub) 2*ones(1,NSub)]; %- timepoint
group = [ones(1,Ngp1) 2*ones(1,Ngp2) ones(1,Ngp1) 2*ones(1,Ngp2)]; %-- group ID for every datapoint

load('beta_closed.mat');
data = beta_closed;
data = mean(data(ele,:),1);

tab = F1LDF1(data,group,subject,time);
results(:,:,1) = tab;


results_pval = results(:,3,:); % '3' for t values - use '2' for n values
p = reshape(results_pval,1,3);
q=0.05;

%%%% FDR correction
%[pID,pN] = fdr(p,q)
%sig = find(p<=pID)

%%%% no FDR correction
sig = find(p<=q)
p(sig)