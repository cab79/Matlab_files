dbstop if error
clear all
close all

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\raincloud_plots');
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
[cb] = cbrewer('qual', 'Set1', 10, 'pchip');

% Load Bayes optimal results
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; 
fname = 'CORE_fittedparameters_percmodel10_bayesopt_20180805T093044.mat';
ls=load(fullfile(S.path.hgf,'fitted',fname));
D_fit=ls.D_fit; 

% traj to plot
S.traj{1} = {
     {'PL'},{'dau','da'},{[0 0]},{[],[]};
     }; 
[S] = HGF_traj2mat(S,D_fit);

% Pred data transformation
S.pred_transform = 'rank'; % arcsinh, rank or notrans
if ~isempty(S.pred_transform)
    for pr = 1:size(S.pred,2)
        if strcmp(S.pred_transform,'arcsinh') % need to modify this to apply to only selected predictors
            x=S.pred(:,pr);
            S.pred(:,pr)=log(x+sqrt(x.^2+1));

        elseif strcmp(S.pred_transform,'rank') % need to modify this to apply to only selected predictors
            x=S.pred(:,pr);
            [~,S.pred(:,pr)]=sort(x);
        end
    end
end
    
% plot distributions as rainclouds
figure
nplot = size(S.pred,2);
for i = 1:nplot
    subplot(nplot,1,i), raincloud_plot(S.pred(:,i), cb(i,:));
    title(S.pred_label{i})
end