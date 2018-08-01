function [output, model] = cosmo_gpml_CAB(samples_train, targets_train, samples_test, opt)

% predicted=cosmo_gpml_CAB(samples_train, targets_train, samples_test[,opt])
%
% Inputs:
%   samples_train       PxR training data for P samples and R features
%   targets_train       Px1 training data classes
%   samples_test        QxR test data
%   opt                 Optional struct with optional field:
%    .xxx               None
%
% Output:
%   predicted          Qx1 predicted data classes for samples_test
%
%-------------------------------------------------------------------------%
% Notes:

% Two target type cause classification to occur; more than two cause
% regression.

% Bayesian multivariate classification or regression using Gaussian Processes
% Requires GPML toolbox as implemented in PRoNTo:
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
% http://www.mlnl.cs.ucl.ac.uk/pronto/

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

% set option defaults
if nargin<4 || isempty(opt)
    opt=struct();
end
if ~isfield(opt, 'max_feature_count'), opt.max_feature_count=50000; end
if ~isfield(opt, 'optimise_theta'), opt.optimise_theta=1; end

% support repeated testing on different data after training every time
% on the same data. This is achieved by caching the training data
% and associated model
persistent cached_targets_train;
persistent cached_samples_train;
persistent cached_opt;
persistent cached_model;

if isequal(cached_targets_train, targets_train) && ...
        isequal(cached_opt, opt) && ...
        isequal(cached_samples_train, samples_train)
    % use cache
    model=cached_model;
else
    % train classifier
    model=train(samples_train, targets_train, opt);

    % store model
    cached_targets_train=targets_train;
    cached_samples_train=samples_train;
    cached_opt=opt;
    cached_model=model;
end

% test classifier
output=test(model, samples_test, samples_train);
%predicted = output.predictions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model=train(samples_train, targets_train, opt)
% train GP

[ntrain, nfeatures]=size(samples_train);
ntrain_=numel(targets_train);

if ntrain_~=ntrain
    error(['size mismatch: samples_train has %d rows, '...
            'targets_train has %d values'],ntrain,ntrain_);
end

if nfeatures>opt.max_feature_count
    % compute size requirements for numbers.
    % tyically this is 8 bytes per number, unless single precision
    % is used.
    w=whos('samples_train');
    size_per_number=w.bytes/numel(samples_train);

    % the covariance matrix is nfeatures x nfeatures
    mem_required=nfeatures^2*size_per_number;
    mem_required_mb=round(mem_required)/1e6; % in megabytes

    error(['A large number of features (%d) was found, '...
            'exceeding the safety limit max_feature_count=%d '...
            'for the %s function.\n'...
            'This limit is imposed because computing the '...
            'covariance matrix will require in the order of '...
            '%.0f MB of memory, and inverting it may take '...
            'a long time.\n'...
            'The safety limit can be changed by setting '...
            'the ''max_feature_count'' option to another '...
            'value, but large values ***may freeze or crash '...
            'the machine***.'],...
            nfeatures,opt.max_feature_count,...
            mfilename(),mem_required_mb);
end

classes=fast_vector_unique(targets_train);
nclasses=numel(classes);

if nclasses==2
    type = 'classification';
elseif nclasses>2
    type = 'regression';
else
    error('incorrect number of target types')
end

% configure covariances
%Kss = d.testcov;
%te_lbs = targets_test; % testing labels

% get parameter defaults
[hyp, maxeval, inffunc, meanfunc, covfunc, likfunc, y, mtr] = gp_params(type,targets_train);

% train
if opt.optimise_theta
    %if map
    %    [hyp,nlmls] = minimize(hyp, @gp_map, maxeval, inffunc, meanfunc, covfunc, likfunc, K, y, priors);
    %else
        [hyp nlmls] = minimize(hyp, @prt_gp, maxeval, inffunc, meanfunc, covfunc, likfunc, {samples_train*samples_train'}, y);
    %end
else
    nlmls = prt_gp(hyp, inffunc, meanfunc, covfunc, likfunc, {samples_train*samples_train'}, y);
end

% Outputs
% -------------------------------------------------------------------------
model=struct();
model.type        = type;
model.loghyper    = hyp;
model.nlml        = min(nlmls);
model.classes=classes;
model.targets_train = targets_train;
model.mtr = mtr;

function output=test(model, samples_test, samples_train)
% make GP predictions
[hyp, maxeval, inffunc, meanfunc, covfunc, likfunc, y, mtr] = gp_params(model.type,model.targets_train);
[ymu,ys2,fmu,fs2,lp,post] = prt_gp(hyp, inffunc, meanfunc, covfunc, likfunc, {samples_train*samples_train'}, y, {samples_test*samples_train'}, zeros(size(samples_test,1),1),{samples_test*samples_test'});

if strcmp(model.type,'classification')
    p = exp(lp);
    output.predictions = (1-real(p > 0.5)) + 1;
    output.func_val    = p;
else % regression
    output.predictions = ymu + model.mtr;
    output.func_val    = output.predictions;
end
output.mu          = ymu;
output.s2          = ys2;
output.alpha       = post.alpha;
%output.sW          = post.sW;
%output.L           = post.L;

% compute weights
nclasses=1; % regression
for c = 1:nclasses
    img1d     = zeros(size(samples_train(1,:)),'single');
    for i=1:length(output.alpha)
        tmp1 = single(samples_train(i,:));
        tmp2 = single(output.alpha(i,c));
        img1d = img1d + tmp1 * tmp2;
    end
    output.class_weight{c} = img1d;
end



function unq_xs=fast_vector_unique(xs)
xs_sorted=sort(xs(:));
idxs=([true;diff(xs_sorted)>0]);
unq_xs=xs_sorted(idxs);

function [hyp, maxeval, inffunc, meanfunc, covfunc, likfunc, y, mtr] = gp_params(type,targets_train)

% Info from PRoNTo prt_gp function:
% Gaussian Process inference and prediction. The gp function provides a
% flexible framework for Bayesian inference and prediction with Gaussian
% processes for scalar targets, i.e. both regression and binary
% classification. The prior is Gaussian process, defined through specification
% of its mean and covariance function. The likelihood function is also
% specified. Both the prior and the likelihood may have hyperparameters
% associated with them.
%
% Two modes are possible: training or prediction: if no test cases are
% supplied, then the negative log marginal likelihood and its partial
% derivatives w.r.t. the hyperparameters is computed; this mode is used to fit
% the hyperparameters. If test cases are given, then the test set predictive
% probabilities are returned. Usage:
%
%   training: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, x, y);
% prediction: [ymu ys2 fmu fs2   ] = gp(hyp, inf, mean, cov, lik, x, y, xs);
%         or: [ymu ys2 fmu fs2 lp] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);
%
% where:
%
%   hyp      column vector of hyperparameters
%   inf      function specifying the inference method 
%   cov      prior covariance function (see below)
%   mean     prior mean function
%   lik      likelihood function
%   x        n by D matrix of training inputs
%   y        column vector of length n of training targets
%   xs       ns by D matrix of test inputs
%   ys       column vector of length nn of test targets
%
%   nlZ      returned value of the negative log marginal likelihood
%   dnlZ     column vector of partial derivatives of the negative
%               log marginal likelihood w.r.t. each hyperparameter
%   ymu      column vector (of length ns) of predictive output means
%   ys2      column vector (of length ns) of predictive output variances
%   fmu      column vector (of length ns) of predictive latent means
%   fs2      column vector (of length ns) of predictive latent variances
%   lp       column vector (of length ns) of log predictive probabilities
%
%   post     struct representation of the (approximate) posterior
%            3rd output in training mode and 6th output in prediction mode


% configure default parameters for GP optimisation (see PRoNTo 2.0 toolbox for these settings)
meanfunc  = @meanConstcell;
covfunc   = @covLINkcell; 
maxeval   = -20;
switch type
    case 'classification'
          %meanfunc = [];                    % empty: don't use a mean function
          %covfunc = @covSEiso;              % Squared Exponental covariance function, see help covSEiso
          likfunc   = @likErf;
          inffunc   = @prt_infEP;
          %hyp = struct('mean', [], 'cov', [0 0], 'lik', -1); % initialize the hyperparameter struct
          mtr       = nan;

    case 'regression'
          %meanfunc = [];                    % empty: don't use a mean function
          %covfunc = @covSEiso;              % Squared Exponental covariance function, see help covSEiso
          likfunc = @likGauss;              % Gaussian likelihood
          inffunc   = @prt_infExact;
          %hyp = struct('mean', [], 'cov', [0 0], 'lik', -1); % initialize the hyperparameter struct
          mtr       = mean(targets_train);      % mean of the training data
end

% check we can reach the binary library
if ~exist('gp','file')
    error(['Error:'...
        ' ''gp'' function could not be found !' ...
        ' SOLUTION: Please check your path.']);
end

% Set default hyperparameters
% -------------------------------------------------------------------------
nhyp = str2num([feval(covfunc); feval(likfunc); feval(meanfunc)]);
if nhyp(1) > 0
    hyp.cov = zeros(nhyp(1),1);
end
if nhyp(2) > 0
    hyp.lik = zeros(nhyp(2),1);
end
if nhyp(3) > 0 
    hyp.mean = zeros(nhyp(3),1);
end

% configure targets
if strcmp(type,'classification')
    % convert targets to +1/-1
    y = -1*(2 * targets_train - 3);
else
    y = targets_train - mtr;
end
    