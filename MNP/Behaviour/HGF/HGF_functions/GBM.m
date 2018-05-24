function [traj, infStates] = GBM(r, p, varargin)
% GENERAL BINARY (INPUT) MODEL

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
    p = GBM_transp(r, p);
end

% Number of levels
l = r.c_prc.n_levels;
try
    %in = r.c_prc.n_inputs; % num input types / conditions
    N = r.c_prc.n_targets; % number of stimulus types that can be a target for each condition. For a binary model, N=1, otherwise a categorical model is needed.
catch
    %in = 1; 
    N=1;
    %l = 100;%(length(p)-3*in-1)/nLevFac;
    
    %if l ~= floor(l)
    %    error('tapas:hgf:UndetNumLevels', 'Cannot determine number of levels');
    %end
end


% Unpack parameters

%CAB: names and indices
nme=r.c_prc.pnames;
idx=r.c_prc.priormusi;

for pn=1:length(nme)
    nme2 = strsplit(nme{pn,1},'log');
    eval([nme2{end} ' = p(idx{pn});']);
end

if exist('om','var')
    if length(om)>1
        th   = exp(om(end));
        om(end)=[];
    end
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
v     = NaN(n,l);
w     = NaN(n,l-1);
da    = NaN(n,l);
dau   = NaN(n,1);
al   = NaN(n,1);

% Representation priors
% Note: first entries of the other quantities remain
% NaN because they are undefined and are thrown away
% at the end; their presence simply leads to consistent
% trial indices.
mu(1,1) = tapas_sgm(mu_0(1), 1);
pi(1,1) = Inf;
if l>1
    mu(1,2:end) = mu_0(2:end);
    pi(1,2:end) = 1./sa_0(2:end);
end
if exist('g_0','var')
    if r.c_prc.one_alpha
        g  = NaN(n,1); % Kalman gain (optional)
    else
        g  = NaN(n,2); % Kalman gain (optional)
    end
    g(1,:)  = g_0; % Kalman gain (optional)
    expom = exp(om);
end

for k=2:1:n
    
    if not(ismember(k-1, r.ign))
        
        if strcmp(r.c_prc.priorsource,'hierarchical')
            if strcmp(r.c_prc.priortype,'fixed')
                if r.c_prc.n_muin>1
                    muhat(1,2) = mu_0(1+u(k,2));
                    pihat(1,2) = 1./sa_0(1+u(k,2));
                else
                    muhat(1,2) = mu_0(2);
                    pihat(1,2) = 1./sa_0(2);
                end
                % 2nd level prediction
                muhat(k,2) = muhat(1,2);%mu(k-1,2) +t(k) *rho(2); % fixed to initial value - not updated on each trial
                
            elseif strcmp(r.c_prc.priortype,'dynamic')
                % 2nd level prediction
                muhat(k,2) = mu(k-1,2) +t(k) *rho(2);
                
            end
            % Prediction from level 2 (which can be either fixed or dynamic)
            muhat(k,1) = tapas_sgm(muhat(k,2), 1);
            
        elseif strcmp(r.c_prc.priorsource,'state')
            % Prediction from prior state, e.g. Kalman filter
            muhat(k,1) =  mu(k-1,1);
            
        end
        
        % Precision of prediction
        pihat(k,1) = 1/(muhat(k,1)*(1 -muhat(k,1)));
        
        
        %% Updates
        
        % Value prediction error, e.g. for Kalman filter, also known as the
        % "innovation"
        dau(k) = u(k,1) -muhat(k,1);
        
        % set alpha
        if r.c_prc.one_alpha
            al(k)=al0(u(k,2));
        else
            if u(k,1)==0
                al(k)=al0(u(k,2));
            elseif u(k,1)==1
                al(k)=al1(u(k,2));
            end
        end
        
        % 
        if strcmp(r.c_prc.priorsource,'hierarchical')
            % Likelihood functions for binary model / SDT: one for each
            % possible signal
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
            mu(k,1) = muhat(k,1) *und1 /(muhat(k,1) *und1 +(1 -muhat(k,1)) *und0);
            
            % Representation prediction error
            da(k,1) = mu(k,1) -muhat(k,1);
            
            % second level predictions and precisions
            if strcmp(r.c_prc.priortype,'fixed')
                mu(k,2) = muhat(k,2); % for a model with higher level predictions, which are fixed
                % At second level, assume Inf precision for a model with invariable predictions
                pi(k,2) = Inf;
                pihat(k,2) = Inf;
                
            elseif strcmp(r.c_prc.priortype,'dynamic')
                % Precision of prediction
                pihat(k,2) = 1/(1/pi(k-1,2) +exp(ka(2) *mu(k-1,3) +om(2)));

                % Updates
                pi(k,2) = pihat(k,2) +1/pihat(k,1);
                mu(k,2) = muhat(k,2) +1/pi(k,2) *da(k,1);

                % Volatility prediction error
                da(k,2) = (1/pi(k,2) +(mu(k,2) -muhat(k,2))^2) *pihat(k,2) -1;
            end
            
            % Implied posterior precision at first level
            sgmmu2 = tapas_sgm(mu(k,2), 1);
            pi(k,1) = pi(k,2)/(sgmmu2*(1-sgmmu2));
            
            if l > 3
                % Pass through higher levels
                % ~~~~~~~~~~~~~~~~~~~~~~~~~~
                for j = 3:l-1
                    % Prediction
                    muhat(k,j) = mu(k-1,j) +t(k) *rho(j);

                    % Precision of prediction
                    pihat(k,j) = 1/(1/pi(k-1,j) +t(k) *exp(ka(j) *mu(k-1,j+1) +om(j)));

                    % Weighting factor
                    v(k,j-1) = t(k) *exp(ka(j-1) *mu(k-1,j) +om(j-1));
                    w(k,j-1) = v(k,j-1) *pihat(k,j-1);

                    % Updates
                    pi(k,j) = pihat(k,j) +1/2 *ka(j-1)^2 *w(k,j-1) *(w(k,j-1) +(2 *w(k,j-1) -1) *da(k,j-1));

                    if pi(k,j) <= 0
                        error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                    end

                    mu(k,j) = muhat(k,j) +1/2 *1/pi(k,j) *ka(j-1) *w(k,j-1) *da(k,j-1);

                    % Volatility prediction error
                    da(k,j) = (1/pi(k,j) +(mu(k,j) -muhat(k,j))^2) *pihat(k,j) -1;
                end
            end
            if l>2
                % Last level
                % ~~~~~~~~~~
                % Prediction
                muhat(k,l) = mu(k-1,l) +t(k) *rho(l);

                % Precision of prediction
                pihat(k,l) = 1/(1/pi(k-1,l) +t(k) *th);

                % Weighting factor
                v(k,l)   = t(k) *th;
                v(k,l-1) = t(k) *exp(ka(l-1) *mu(k-1,l) +om(l-1));
                w(k,l-1) = v(k,l-1) *pihat(k,l-1);

                % Updates
                pi(k,l) = pihat(k,l) +1/2 *ka(l-1)^2 *w(k,l-1) *(w(k,l-1) +(2 *w(k,l-1) -1) *da(k,l-1));

                if pi(k,l) <= 0
                    error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                end

                mu(k,l) = muhat(k,l) +1/2 *1/pi(k,l) *ka(l-1) *w(k,l-1) *da(k,l-1);

                % Volatility prediction error
                da(k,l) = (1/pi(k,l) +(mu(k,l) -muhat(k,l))^2) *pihat(k,l) -1;
            end

        elseif strcmp(r.c_prc.priorsource,'state') % Kalman
            % Gain update - optimal gain is calculated from ratio of input
            % variance to representation variance
            
            % Same gain function modified by two different alphas
            g(k) = (g(k-1) +al(k)*expom)/(g(k-1) +al(k)*expom +1);
            % Hidden state mean update
            mu(k,1) = muhat(k,1)+g(k)*dau(k);
            pi(k,1) = (1-g(k)) *al(k)*expom; 
            
            % Alternative: separate gain functions for each stimulus type
      %      if r.c_prc.one_alpha
      %          pi_u=al0(u(k,2));
      %          g(k,1) = (g(k-1,1) +pi_u*expom)/(g(k-1,1) +pi_u*expom +1);
      %          % Hidden state mean update
      %          mu(k,1) = muhat(k,1)+g(k,1)*dau(k);
      %          pi(k,1) = (1-g(k,1)) *pi_u*expom;
      %      else
      %          if u(k,1)==0
      %              pi_u=al0(u(k,2));
      %              g(k,1) = (g(k-1,1) +pi_u*expom)/(g(k-1,1) +pi_u*expom +1);
      %              g(k,2) = g(k-1,2);
      %              % Hidden state mean update
      %              mu(k,1) = muhat(k,1)+g(k,1)*dau(k);
      %              pi(k,1) = (1-g(k,1)) *pi_u*expom;
      %          elseif u(k,1)==1
      %              pi_u=al1(u(k,2));
      %              g(k,2) = (g(k-1,2) +pi_u*expom)/(g(k-1,2) +pi_u*expom +1);
      %              g(k,1) = g(k-1,1);
      %              % Hidden state mean update
      %              mu(k,1) = muhat(k,1)+g(k,2)*dau(k);
      %              pi(k,1) = (1-g(k,2)) *pi_u*expom;
      %          end
      %      end
      
            % Representation prediction error
            da(k,1) = mu(k,1) -muhat(k,1);
            
        end
        
        % RESPONSE BIAS
        if exist('rb','var')
            mu(k,1) = mu(k,1)+rb;
        end
        
        
        
    else
        mu(k,:) = mu(k-1,:); 
        pi(k,:) = pi(k-1,:);

        muhat(k,:) = muhat(k-1,:);
        pihat(k,:) = pihat(k-1,:);
        
        v(k,:)  = v(k-1,:);
        w(k,:)  = w(k-1,:);
        da(k,:) = da(k-1,:);
        dau(k) = dau(k-1);
        if exist('g','var')
            g(k,:)=g(k-1,:);
        end;
        al(k)  = al(k-1);
    end
end

if l>1
    % Implied learning rate at the first level
    sgmmu2 = tapas_sgm(mu(:,2), 1);
    lr1    = diff(sgmmu2)./da(2:n,1);
    lr1(da(2:n,1)==0) = 0;
end

% Remove representation priors
mu(1,:)  = [];
pi(1,:)  = [];
% Check validity of trajectories
if any(isnan(mu(:))) %|| any(isnan(pi(:)))
    error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
else
    % Check for implausible jumps in trajectories
    % CAB: only use first 500 trials - after that changes in precision become too small
    ntrials = min(length(mu),500);
    dmu = diff(mu(1:ntrials,2:end));
    dpi = diff(pi(1:ntrials,2:end));
    rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
    rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

    jumpTol = 16;
    if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
        disp('Use plot for diagnosis: see within function'); % plot(abs(dpi(:,2))); hold on; plot(rmdpi(:,2),'r'); hold on; plot(jumpTol*rmdpi(:,2),'g')
    end
end

% Remove other dummy initial values
muhat(1,:) = [];
pihat(1,:) = [];
v(1,:)     = [];
w(1,:)     = [];
da(1,:)    = [];
dau(1)     = [];
if exist('g','var')
    g(1,:)  = [];
end
al(1)     = [];

% Create result data structure
traj = struct;

traj.mu     = mu;
traj.sa     = 1./pi;

traj.muhat  = muhat;
traj.sahat  = 1./pihat;
traj.v      = v;
traj.w      = w;
traj.da     = da;
traj.dau    = dau;
if exist('g','var')
    traj.g     = g;
end

% Updates with respect to prediction
traj.ud = muhat -mu;

% Psi (precision weights on prediction errors)
psi        = NaN(n-1,l+1);
psi(:,1)   = 1./(al.*pi(:,1)); % dot multiply only if al is a vector
if l>1; psi(:,2)   = 1./pi(:,2);end
if l>2; psi(:,3:l) = pihat(:,2:l-1)./pi(:,3:l);end
traj.psi   = psi;

% Epsilons (precision-weighted prediction errors)
epsi        = NaN(n-1,l);
epsi(:,1)   = psi(:,1) .*dau;
if l>1; epsi(:,2:l) = psi(:,2:l) .*da(:,1:l-1);end
traj.epsi   = epsi;

% Full learning rate (full weights on prediction errors)
wt        = NaN(n-1,l);
if l==1; lr1=psi(:,1);end
wt(:,1)   = lr1;
if l>1; wt(:,2)   = psi(:,2); end
if l>2; wt(:,3:l) = 1/2 *(v(:,2:l-1) *diag(ka(2:l-1))) .*psi(:,3:l); end
traj.wt   = wt;

% Create matrices for use by the observation model
infStates = NaN(n-1,l,7);
infStates(:,:,1) = traj.muhat;
infStates(:,:,2) = traj.sahat;
infStates(:,:,3) = traj.mu;
infStates(:,:,4) = traj.sa;
infStates(:,:,5) = traj.da;
infStates(:,:,6) = traj.epsi;
infStates(:,1,7) = traj.dau;

return;
