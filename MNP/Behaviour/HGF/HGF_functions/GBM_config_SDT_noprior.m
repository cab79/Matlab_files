function c = GBM_config_SDT_noprior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Config structure
c = struct;

% Model name
c.model = 'GBM';

% Number of levels
c.n_levels = 2;

% Number of inputs (minimum: 1)
c.n_alphain = 2;
c.n_muin = 1;

% Number of target stimuli
c.n_targets = 1;

% Source of prior
c.priorsource='hierarchical'; % from one level above in the hierarchy. Higher level can have a fixed or dynamic prior.
%c.priorsource='state'; % from the previous state at time t-1. e.g. Kalman Filter.

% Prior type
c.priortype='fixed';
%c.priortype='dynamic';

% Input intervals
% If input intervals are irregular, the last column of the input
% matrix u has to contain the interval between inputs k-1 and k
% in the k-th row, and this flag has to be set to true
c.irregular_intervals = false;

% Sufficient statistics of Gaussian parameter priors

% Initial mus and sigmas
% Format: row vectors of length n_levels

%% Initial hidden state mean and variance

% HIERARCHICAL MODELS
% For all but the first two levels, this is usually best
% kept fixed to 1 (determines origin on x_i-scale). The 
% first level is NaN because it is determined by the second,
% and the second implies neutrality between outcomes when it
% is centered at 0.

% Dynamic Priors
%c.mu_0mu = [NaN, 0, 1];
%c.mu_0sa = [NaN, 0, 0];
%c.logsa_0mu = [NaN,   log(0.1), log(1)];
%c.logsa_0sa = [NaN,          0,      0];

% Fixed Priors
c.mu_0mu = [NaN, repmat(0.5,1,c.n_muin)];
c.mu_0sa = [NaN, repmat(0,1,c.n_muin)];% fixed prior
%c.mu_0sa = [NaN, repmat(1,1,c.n_muin)];% estimated prior
c.logsa_0mu = [NaN, repmat(log(0.1),1,c.n_muin)];
c.logsa_0sa = [NaN, repmat(0,1,c.n_muin)];

% PRIOR STATE MODELS
% e.g. Kalman Filter or Rascorla-Wagner
% 99991   Value of the first input
%         Usually a good choice for mu_0mu(1)
% 99992   Variance of the first 20 inputs
%         Usually a good choice for mu_0sa(1)
%c.mu_0mu = 99991;
%c.mu_0sa = 99992;

%% Input variance
% Alpha
% Format: scalar.
c.logal0mu = repmat(log(0.5),1,c.n_alphain);
c.logal0sa = repmat(1,1,c.n_alphain);
c.logal1mu = repmat(log(0.5),1,c.n_alphain);
c.logal1sa = repmat(0,1,c.n_alphain); % fixed
c.one_alpha =1; % 1 = two alphas are the same. Must fixed variance of al1 prior to reduce complexity.

%% Gain

% STATE MODELS
% Initial Kalman gain
c.logg_0mu = 0.1;
c.logg_0sa = 1;

% HIERARCHICAL MODELS

% Rhos
% Format: row vector of length n_levels.
% Undefined (therefore NaN) at the first level.
% Fix this to zero to turn off drift.
c.rhomu = [NaN, 0, 0];
c.rhosa = [NaN, 0, 0];

% Kappas
% Format: row vector of length n_levels-1.
% Undefined (therefore NaN) at the first level.
% This should be fixed (preferably to 1) if the observation model
% does not use mu_i+1 (kappa then determines the scaling of x_i+1).
c.logkamu = [NaN, log(1)];
c.logkasa = [NaN,      0];

%% Process variance

% HIERARCHICAL MODELS
% Format: row vector of length n_levels.
% Undefined (therefore NaN) at the first level.
c.ommu = [NaN,  -6,  -6];
c.omsa = [NaN, 4^2, 4^2];

% STATE MODELS (for Kalman Filter)
% 99993   Log-variance of the first 20 inputs
%         Usually a good choice for logsa_0mu(1), and
%         its negative, ie the log-precision of the
%         first 20 inputs, for logpiumu
%c.ommu = 99993;
%c.omsa = 1;

%% Expected input means
% HIERARCHICAL MODELS ONLY
% Eta0
% Format: scalar.
c.eta0mu = 0;
c.eta0sa = 0;

% Eta1
% Format: scalar.
c.eta1mu = 1;
c.eta1sa = 0;

%%
% Response bias
% Format: scalar.
c.rbmu = 0;
c.rbsa = 0;

%% Parameter names
% only names added here will be estimated
% second column: 1 if it is a variance parameter
% SDT MODELS
c.pnames = {
    'mu_0',0
    'logsa_0',1
    'logal0',1
    'logal1',1
    'eta0',0
    'eta1',0
    'rb',0
    };

% KF MODELS
%c.pnames = {
%    'mu_0',0
%    'logal0',1
%    'logal1',1
%    'rb',0
%    'logg_0',1
%    'om',0
%    };

% kf MODELS
%c.pnames = {
%    'mu_0',0
%    'logsa_0',1
%    'logal0',1
%    'logal1',1
%    'rho',0
%    'logka',1
%    'om',0
%    'eta0',0
%    'eta1',0
%    };

%% Gather prior settings and their indices in vectors
nparams =[];
c.priormus=[];
c.priorsas=[];
for pn = 1:length(c.pnames)
    eval(['c.priormus = [c.priormus c.' c.pnames{pn,1} 'mu];']);
    eval(['c.priorsas = [c.priorsas c.' c.pnames{pn,1} 'sa];']);
    eval(['nparams(pn) = length(c.' c.pnames{pn,1} 'mu);']);
    if pn==1
        st = 0;
    else
        st=sum(nparams(1:pn-1));
    end
    c.priormusi{pn} = st+1:sum(nparams(1:pn));
end

% Check whether we have the right number of priors
%expectedLength = 2*c.n_levels+c.n_inputs*4+1;
%if length([c.priormus, c.priorsas]) ~= 2*expectedLength;
%    error('tapas:hgf:PriorDefNotMatchingLevels', 'Prior definition does not match number of levels.')
%end

%% Model function handle
c.prc_fun = @GBM;

% Handle to function that transforms perceptual parameters to their native space
% from the space they are estimated in
c.transp_prc_fun = @GBM_transp;

return;
