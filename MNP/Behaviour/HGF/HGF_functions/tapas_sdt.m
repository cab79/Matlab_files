function [traj, infStates] = tapas_sdt(r, p, varargin)

% based on "evaluateYesNo" model from https://link.springer.com/article/10.3758/s13414-014-0830-0

% NOTES
% 1. in and N defined separately
% 2. representation Mu has no variance in SDT, there is just noise (alpha)
% that influences it's mean value.
% 3. must be a two-level model because the 2nd level sets the prior for the
% first, even if the prior is fixed

% PLANNED UPDATES
% 1. Allow separate evaluation of alpha for targets and non-targets - done
% 2. Estimate prior expectation (fixed)
% 3. Prior and alpha variable by block type
% 4. Prior and alpha variable by trial (Kalman)
% 5. Hierarchical priors

% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans');
    p = tapas_sdt_transp(r, p);
end

nLevFac = 2; % number of parameters that have length "r.c_prc.n_levels"
%nInpFac = 3; % number of parameters that have length "r.c_prc.n_inputs"

% Number of levels
try
    l = r.c_prc.n_levels;
    %in = r.c_prc.n_inputs; % num input types / conditions
    N = r.c_prc.n_targets; % number of stimulus types that can be a target for each condition. For a binary model, N=1, otherwise a categorical model is needed.
catch
    %in = 1; 
    N=1;
    l = 100;%(length(p)-3*in-1)/nLevFac;
    
    if l ~= floor(l)
        error('tapas:hgf:UndetNumLevels', 'Cannot determine number of levels');
    end
end


% Unpack parameters

%CAB: names and indices
nme=r.c_prc.pnames;
idx=r.c_prc.priormusi;

for pn=1:length(nme)
    nme2 = strsplit(nme{pn,1},'log');
    eval([nme2{end} ' = p(idx{pn});']);
end

%mu_0 = p(1:l);
%sa_0 = p(l+1:2*l);
%al0   = p(nLevFac*l+1:nLevFac*l+in);
%al1   = p(nLevFac*l+in+1:nLevFac*l+2*in);
%eta0 = p(nLevFac*l+2*in+1:nLevFac*l+3*in);
%eta1 = p(nLevFac*l+3*in+1:nLevFac*l+4*in);
%rb = p(nLevFac*l+4*in+1);

% Add dummy "zeroth" trial
u = [zeros(1,2); r.u(:,1:2)];

% Number of trials (including prior)
n = size(u,1);

% Assume that if u has more than one column, the last contains t
try
    if r.c_prc.irregular_intervals
        if size(u,2) > 1
            t = [0; r.u(:,end)];
        else
            error('tapas:hgf:InputSingleColumn', 'Input matrix must contain more than one column if irregular_intervals is set to true.');
        end
    else
        t = ones(n,1);
    end
catch
    if size(u,2) > 1
        t = [0; r.u(:,end)];
    else
        t = ones(n,1);
    end
end

% Initialize updated quantities

% Representations
mu = NaN(n,l);
pi = NaN(n,l);

% Other quantities
muhat = NaN(n,l);
pihat = NaN(n,l);

% Representation priors
% Note: first entries of the other quantities remain
% NaN because they are undefined and are thrown away
% at the end; their presence simply leads to consistent
% trial indices.
mu(1,1) = tapas_sgm(mu_0(1), 1);
pi(1,1) = Inf;

for k=2:1:n
    
    if not(ismember(k-1, r.ign))
        
        if r.c_prc.n_muin>1
            muhat(1,2) = mu_0(1+u(k,2));
            pihat(1,2) = 1./sa_0(1+u(k,2));
        else
            muhat(1,2) = mu_0(2);
            pihat(1,2) = 1./sa_0(2);
        end
	
        % 2nd level prediction
        muhat(k,2) = muhat(1,2);%mu(k-1,2) +t(k) *rho(2); % fixed to initial value - not updated on each trial

        % 1st level
        % ~~~~~~~~~
        % Prediction
        muhat(k,1) = tapas_sgm(muhat(k,2), 1);
        
        % Precision of prediction
        pihat(k,1) = 1/(muhat(k,1)*(1 -muhat(k,1)));
        
        % Mean update
        mu(k,1) = u(k,1);
        
        % after equating mu to u, it does not need updating unless u
        % differs from eta (mean of true signal)
        if r.c_prc.n_alphain >1
            if r.c_prc.one_alpha
                und1 = exp(-(u(k,1) -eta1)^2/(2*al0(u(k,2))));
            else
                und1 = exp(-(u(k,1) -eta1)^2/(2*al1(u(k,2))));
            end
            und0 = exp(-(u(k,1) -eta0)^2/(2*al0(u(k,2))));
        else
            if r.c_prc.one_alpha
                und1 = exp(-(u(k,1) -eta1)^2/(2*al0));
            else
                und1 = exp(-(u(k,1) -eta1)^2/(2*al1));
            end
            und0 = exp(-(u(k,1) -eta0)^2/(2*al0));
        end
        %und1 = exp(-(u(k) -eta1(1))^2/(2*al(1)));
        %und0 = exp(-(u(k) -eta0(1))^2/(2*al(1)));
        %eta0(:)=eta0(1);
        %eta1(:)=eta1(1);
        %al(:)=al(1);
        
        %% Updates
        
        mu(k,1) = muhat(k,1) *und1 /(muhat(k,1) *und1 +(1 -muhat(k,1)) *und0);
        
        %% ALTERNATIVE

        % likelihood of D
        % normal distribution describing likelihood of the input
    %    Ld0 = prod( normpdf(u(k,1), eta0(u(k,2)), al(u(k,2))) );
    %    Ld1 = prod( normpdf(u(k,1), eta1(u(k,2)), al(u(k,2))) );%

    %    PosteriorD0 = Ld0 .* (1-muhat(k,1));					% posterior
    %    PosteriorD1 = Ld1 .* muhat(k,1);					% posterior
    %    mu(k,1) = PosteriorD1./(PosteriorD0+PosteriorD1);	% normalise
        
        
        
        %% STEP 3: RESPONSE BIAS

        mu(k,1) = mu(k,1)+rb;
        
        % Equate 2nd level representation to the prior, as it doesn't
        % change
        mu(k,2) = muhat(k,1);
        
        % Assume Inf precision
        pi(k,:) = Inf;
        pihat(k,2) = Inf;

    else
        mu(k,:) = mu(k-1,:); 
        pi(k,:) = pi(k-1,:);

        muhat(k,:) = muhat(k-1,:);
        pihat(k,:) = pihat(k-1,:);
        
    end
	
end


% Remove representation priors
mu(1,:)  = [];
pi(1,:)  = [];

% Check validity of trajectories
if any(isnan(mu(:))) %|| any(isnan(pi(:)))
    error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
else
    % Check for implausible jumps in trajectories
    dmu = diff(mu(:,2:end));
    dpi = diff(pi(:,2:end));
    rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
    rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

    jumpTol = 16;
    if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
    end
end

% Remove other dummy initial values
muhat(1,:) = [];
pihat(1,:) = [];

% Create result data structure
traj = struct;

traj.mu     = mu;
traj.sa     = 1./pi;

traj.muhat  = muhat;
traj.sahat  = 1./pihat;

% Updates with respect to prediction
traj.ud = muhat -mu;

% Psi (precision weights on prediction errors)
psi        = NaN(n-1,l);
psi(:,2)   = 1./pi(:,2);

% Epsilons (precision-weighted prediction errors)
epsi        = NaN(n-1,l);
traj.epsi   = epsi;

% Create matrices for use by the observation model
infStates = NaN(n-1,l,4);
infStates(:,:,1) = traj.muhat;
infStates(:,:,2) = traj.sahat;
infStates(:,:,3) = traj.mu;
infStates(:,:,4) = traj.sa;

return;
