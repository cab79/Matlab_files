addpath(genpath('C:\Data\Matlab\spm12'));

% data and predictors: all normalized to zero-mean and unit-variance (z score).
n=500; % number of data points
p=1; % number of parameters
y = [1:500]'; % dependent variable (to predict)

%enforces positively constraints on the covariance hyperparameters
% by adopting a log-normal [flat] hyperprior.
lognormhyper=0;

% To implement non-hierarchical Bayes with priors on the parameters use
% a two level model setting the second level design matrix to zeros, 
% providing an unconstrained prior on the first-level parameters
% which allows for single-level Bayesian inference
P{1}.X = [1:500]'; %as your [n x p] design matrix
P{1}.C{1}=eye(n); %as your covariance structure on observation errors
P{2}.X=zeros(p,1);
P{2}.C{2}=eye(p); %as your prior covariance on parameters
[C,P,F] = spm_PEB(y,P,lognormhyper)
scatter(P{1}.X,y)

% For studies with L different subjects  ? {1, ..., L}, the second level can be used
% to facilitate a mixed effects analysis. 
% 1st level: participant-specific parameters. 
% 2nd level: as deviations from the corresponding group parameters. 
% Each subject is modeled with an individual set of first-level parameters ?(1) which are shrunk
% towards subject-independent second level group parameters ?(2)
% Then, use an unconstrained prior on the second level parameters ?(2) by employing an
% all-zero third-level design matrix. The third level functions as a
% shrinkage prior on the group parameters.