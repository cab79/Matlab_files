function c = tapas_sdt_config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Config structure
c = struct;

% Model name
c.model = 'tapas_sdt';

% Number of levels (minimum: 2)
c.n_levels = 2;
% Number of inputs (minimum: 1)
c.n_alphain = 4;
c.n_muin = 4;
% Number of target stimuli
c.n_targets = 1;

% Input intervals
% If input intervals are irregular, the last column of the input
% matrix u has to contain the interval between inputs k-1 and k
% in the k-th row, and this flag has to be set to true
c.irregular_intervals = false;

% Sufficient statistics of Gaussian parameter priors

% Initial mus and sigmas
% Format: row vectors of length n_levels
% For all but the first two levels, this is usually best
% kept fixed to 1 (determines origin on x_i-scale). The 
% first level is NaN because it is determined by the second,
% and the second implies neutrality between outcomes when it
% is centered at 0.
c.mu_0mu = [NaN, repmat(0,1,c.n_muin)];
%c.mu_0sa = [NaN, repmat(0,1,c.n_inputs)];% fixed prior
c.mu_0sa = [NaN, repmat(1,1,c.n_muin)];% variable prior

c.logsa_0mu = [NaN, repmat(log(0.1),1,c.n_muin)];
c.logsa_0sa = [NaN, repmat(0,1,c.n_muin)];

% Alpha
% Format: scalar.
c.logal0mu = repmat(log(0.5),1,c.n_alphain);
c.logal0sa = repmat(1,1,c.n_alphain);
c.logal1mu = repmat(log(0.5),1,c.n_alphain);
c.logal1sa = repmat(0,1,c.n_alphain); % fixed
c.one_alpha =1; % 1 = two alphas are the same. Must fixed variance of al1 prior to reduce complexity.

% Eta0
% Format: scalar.
c.eta0mu = 0;
c.eta0sa = 0;

% Eta1
% Format: scalar.
c.eta1mu = 1;
c.eta1sa = 0;

% Response bias
% Format: scalar.
c.rbmu = 0;
c.rbsa = 0;

% CAB parameter names
%second column: 1 if it is a variance parameter
c.pnames = {
    'mu_0',0
    'logsa_0',1
    'logal0',1
    'logal1',1
    'eta0',0
    'eta1',0
    'rb',0
    };

% CAB Gather prior settings and their indices in vectors
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

% Model function handle
c.prc_fun = @tapas_sdt;

% Handle to function that transforms perceptual parameters to their native space
% from the space they are estimated in
c.transp_prc_fun = @tapas_sdt_transp;

return;
