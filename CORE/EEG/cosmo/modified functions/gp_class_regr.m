function out = gp_class_regr(dat,targets,type,nfold,ndec)
% Bayesian multivariate classification or regression using Gaussian Processes
% http://www.gaussianprocess.org/gpml/code/matlab/doc/

% A Gaussian Process is fully specified by a mean function and a covariance function. 
% These functions are specified separately, and consist of a specification
% of a functional form as well as a set of parameters called hyperparameters.

% For functional forms try:
% help meanFunctions
% help covFunctions

% The likelihood function specifies the probability of the observations given the GP (and the hyperparameters). 
% Several likelihood functions are available, an overview is provided by the likFunctions help function (type help likFunctions to get help). 
% Example likelihood functions include likGauss the Gaussian likelihood for regression and likLogistic the logistic likelihood for classification.

% Several inference methods are avaiable, and overview is provided by the infMethods help file (type help infMethods to get help). 
% An example inference method is infGaussLik for exact inference (regression with Gaussian likelihood).

dat=round(dat,ndec);
rm=all(isnan(dat),1) | all(diff(dat)==0); % remove nans and constants
dat(:,rm)=[];

if isempty(dat)
    out.empty=1;
    return
end

% classification
switch type
    case 'class'
        
    case 'regr'
          meanfunc = [];                    % empty: don't use a mean function
          covfunc = @covSEiso;              % Squared Exponental covariance function, see help covSEiso
          likfunc = @likGauss;              % Gaussian likelihood
          hyp = struct('mean', [], 'cov', [0 0], 'lik', -1); % initialize the hyperparameter struct
end

% GP usage:
%     training: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, dat, targets);
%   prediction: [ymu ys2 fmu fs2   ] = gp(hyp, inf, mean, cov, lik, dat, targets, xs);
%           or: [ymu ys2 fmu fs2 lp] = gp(hyp, inf, mean, cov, lik, dat, targets, xs, ys);
% dat: training input data
% targets: training output targets: can be either binary (classification) or continuous (regression)
% nlZ is the negative log marginal likelihood and dnlZ its partial derivatives wrt the hyperparameters 
% (which are used for training the hyperparameters). The prediction outputs are ymu and ys2 for test 
% output mean and covariance, and fmu and fs2 are the equivalent quenteties for the corresponding latent variables. 
% Finally, lp are the test output log probabilities.
[nlZ,dnlZ] = gp(hyp, inf, mean, cov, lik, dat, targets);

