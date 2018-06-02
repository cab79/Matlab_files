function [traj] = GBM(r, pvec, varargin)
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
mu=[];
% Transform paramaters back to their native space if needed
if ~isempty(varargin) 
    if strcmp(varargin{1},'trans')
        pvec = GBM_transp(r, pvec);
    end
end

% Add dummy "zeroth" trial
u = [zeros(1,size(r.u,2)); r.u(:,:)];

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

%% INITIALISE

% Create param struct
nme=r.c_prc.pnames;
idx=r.c_prc.priormusi;
type='like';
for pn=1:length(nme)
    if strcmp(nme{pn,1}(1:length(type)),type)
        nme2 = strsplit(nme{pn,1},[type '_']);
        nme3 = strsplit(nme2{end},'log');
        M.(type).p.(nme3{end}) = pvec(idx{pn});
    end
end

% joint trajectories (if needed)
M.like.tr.vj_mu = NaN(n,1);
M.like.tr.xc = NaN(n,1);
M.like.tr.xchat = NaN(n,1);
    
for m=1:r.c_prc.nModels
    type = r.c_prc.type{m};
    
    % Create param struct
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            nme2 = strsplit(nme{pn,1},[type '_']);
            nme3 = strsplit(nme2{end},'log');
            M.(type).p.(nme3{end}) = pvec(idx{pn});
        end
    end

    % Number of levels
    M.(type).l = r.c_prc.(type).n_priorlevels;

    if isfield(M.(type).p,'om')
        if length(M.(type).p.om)>1
            M.(type).p.th   = exp(M.(type).p.om(end));
            M.(type).p.om(end)=[];
        end
    end

    % Representations
    M.(type).tr.mu0 = NaN(n,1);
    M.(type).tr.mu = NaN(n,M.(type).l);
    M.(type).tr.pi = NaN(n,M.(type).l);

    % Other quantities
    M.(type).tr.muhat = NaN(n,M.(type).l);
    M.(type).tr.pihat = NaN(n,M.(type).l);
    M.(type).tr.v     = NaN(n,M.(type).l);
    M.(type).tr.w     = NaN(n,M.(type).l-1);
    M.(type).tr.da    = NaN(n,M.(type).l);
    M.(type).tr.dau   = NaN(n,1);
    M.(type).tr.al   = NaN(n,1);

    % Representation priors
    % Note: first entries of the other quantities remain
    % NaN because they are undefined and are thrown away
    % at the end; their presence simply leads to consistent
    % trial indices.
    M.(type).tr.mu(1,1) = tapas_sgm(M.(type).p.mu_0(1), 1);
    M.(type).tr.pi(1,1) = Inf;
    if M.(type).l>1
        M.(type).tr.mu(1,2:end) = M.(type).p.mu_0(2:end);
        M.(type).tr.pi(1,2:end) = 1./M.(type).p.sa_0(2:end);
    end
    if isfield(M.(type).p,'g_0')
        if r.c_prc.(type).one_alpha
            g  = NaN(n,1); % Kalman gain (optional)
        else
            g  = NaN(n,2); % Kalman gain (optional)
        end
        g(1,:)  = M.(type).p.g_0; % Kalman gain (optional)
        M.(type).p.expom = exp(M.(type).p.om);
    end
end

%% UPDATES
for k=2:1:n
    
    for m=1:r.c_prc.nModels
        if not(ismember(k-1, r.ign))
            
            % Unpack likelihood parameters
            type='like';
            pnames = fieldnames(M.(type).p);
            for pn=1:length(pnames)
                p.(pnames{pn}) = M.(type).p.(pnames{pn});
            end

            type = r.c_prc.type{m};
            l=M.(type).l;
            
            % Unpack prior parameters
            pnames = fieldnames(M.(type).p);
            for pn=1:length(pnames)
                p.(pnames{pn}) = M.(type).p.(pnames{pn});
            end
            
            %Unpack prior traj
            tnames = fieldnames(M.(type).tr);
            for tn=1:length(tnames)
                tr.(tnames{tn}) = M.(type).tr.(tnames{tn});
            end
        
            %% Predictions (from previous trial or fixed parameters)
            if strcmp(r.c_prc.(type).priortype,'hierarchical')
                if strcmp(r.c_prc.(type).priorupdate,'fixed')
                    if r.c_prc.(type).n_muin>1
                        tr.muhat(1,2) = p.mu_0(1+u(k,2));
                        tr.pihat(1,2) = 1./p.sa_0(1+u(k,2));
                    else
                        tr.muhat(1,2) = p.mu_0(2);
                        tr.pihat(1,2) = 1./p.sa_0(2);
                    end
                    % 2nd level prediction
                    tr.muhat(k,2) = tr.muhat(1,2);%mu(k-1,2) +t(k) *rho(2); % fixed to initial value - not updated on each trial

                elseif strcmp(r.c_prc.(type).priorupdate,'dynamic')
                    % 2nd level prediction
                    tr.muhat(k,2) = tr.mu(k-1,2) +t(k) *p.rho(2);

                end
                % Prediction from level 2 (which can be either fixed or dynamic)
                tr.muhat(k,1) = tapas_sgm(tr.muhat(k,2), 1);

            elseif strcmp(r.c_prc.(type).priortype,'state')
                % Prediction from prior state, e.g. Kalman filter
                tr.muhat(k,1) =  tr.mu(k-1,1);

            end

            % Precision of prediction
            tr.pihat(k,1) = 1/(tr.muhat(k,1)*(1 -tr.muhat(k,1)));


            %% Updates

            % Value prediction error, e.g. for Kalman filter, also known as the
            % "innovation"
            tr.dau(k) = u(k,1) -tr.muhat(k,1);

            % set alpha
            if ~isfield(p,'al1')
                p.al1=p.al0;
            end
            if u(k,1)==0
                tr.al(k)=p.al0(u(k,2));
            elseif u(k,1)==1
                tr.al(k)=p.al1(u(k,2));
            end

            % 
            if strcmp(r.c_prc.(type).priortype,'hierarchical')
                % Likelihood functions: one for each
                % possible signal
                if r.c_prc.n_inputcond >1
                    und1 = exp(-(u(k,1) -p.eta1)^2/(2*p.al1(u(k,2))));
                    und0 = exp(-(u(k,1) -p.eta0)^2/(2*p.al0(u(k,2))));
                else
                    und1 = exp(-(u(k,1) -p.eta1)^2/(2*p.al1));
                    und0 = exp(-(u(k,1) -p.eta0)^2/(2*p.al0));
                end

                % ORIGINAL
                %mu(k,1) = muhat(k,1) *und1 /(muhat(k,1) *und1 +(1 -muhat(k,1)) *und0);

                %% MOD: need to update this properly
                if u(k,3)==2
                    tr.mu0(k,1) = tr.muhat(k,1) *und1 /(tr.muhat(k,1) *und1 +(1 -tr.muhat(k,1)) *und0);
                    tr.mu(k,1) = tr.mu0(k,1);
                elseif u(k,3)==1
                    tr.mu0(k,1) = (1-tr.muhat(k,1)) *und1 /(tr.muhat(k,1) *und0 +(1 -tr.muhat(k,1)) *und1);
                    tr.mu(k,1) = 1-tr.mu0(k,1);

                    % calculate prediction error for mu0 - muhat

                end


                %%
                % Representation prediction error
                tr.da(k,1) = tr.mu(k,1) -tr.muhat(k,1);

                % second level predictions and precisions
                if strcmp(r.c_prc.(type).priorupdate,'fixed')
                    tr.mu(k,2) = tr.muhat(k,2); % for a model with higher level predictions, which are fixed
                    % At second level, assume Inf precision for a model with invariable predictions
                    tr.pi(k,2) = Inf;
                    tr.pihat(k,2) = Inf;

                elseif strcmp(r.c_prc.(type).priorupdate,'dynamic')
                    % Precision of prediction
                    tr.pihat(k,2) = 1/(1/tr.pi(k-1,2) +exp(p.ka(2) *tr.mu(k-1,3) +p.om(2)));

                    % Updates
                    tr.pi(k,2) = tr.pihat(k,2) +1/tr.pihat(k,1);
                    tr.mu(k,2) = tr.muhat(k,2) +1/tr.pi(k,2) *tr.da(k,1);

                    % Volatility prediction error
                    tr.da(k,2) = (1/tr.pi(k,2) +(tr.mu(k,2) -tr.muhat(k,2))^2) *tr.pihat(k,2) -1;
                end

                % Implied posterior precision at first level
                sgmmu2 = tapas_sgm(tr.mu(k,2), 1);
                tr.pi(k,1) = tr.pi(k,2)/(sgmmu2*(1-sgmmu2));

                if l > 3
                    % Pass through higher levels
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for j = 3:l-1
                        % Prediction
                        tr.muhat(k,j) = tr.mu(k-1,j) +t(k) *p.rho(j);

                        % Precision of prediction
                        tr.pihat(k,j) = 1/(1/tr.pi(k-1,j) +t(k) *exp(p.ka(j) *tr.mu(k-1,j+1) +p.om(j)));

                        % Weighting factor
                        tr.v(k,j-1) = t(k) *exp(p.ka(j-1) *tr.mu(k-1,j) +p.om(j-1));
                        tr.w(k,j-1) = tr.v(k,j-1) *tr.pihat(k,j-1);

                        % Updates
                        tr.pi(k,j) = tr.pihat(k,j) +1/2 *p.ka(j-1)^2 *tr.w(k,j-1) *(tr.w(k,j-1) +(2 *tr.w(k,j-1) -1) *tr.da(k,j-1));

                        if tr.pi(k,j) <= 0
                            error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                        end

                        tr.mu(k,j) = tr.muhat(k,j) +1/2 *1/tr.pi(k,j) *p.ka(j-1) *tr.w(k,j-1) *tr.da(k,j-1);

                        % Volatility prediction error
                        tr.da(k,j) = (1/tr.pi(k,j) +(tr.mu(k,j) -tr.muhat(k,j))^2) *tr.pihat(k,j) -1;
                    end
                end
                if l>2
                    % Last level
                    % ~~~~~~~~~~
                    % Prediction
                    tr.muhat(k,l) = tr.mu(k-1,l) +t(k) *p.rho(l);

                    % Precision of prediction
                    tr.pihat(k,l) = 1/(1/tr.pi(k-1,l) +t(k) *p.th);

                    % Weighting factor
                    tr.v(k,l)   = t(k) *p.th;
                    tr.v(k,l-1) = t(k) *exp(p.ka(l-1) *tr.mu(k-1,l) +p.om(l-1));
                    tr.w(k,l-1) = tr.v(k,l-1) *tr.pihat(k,l-1);

                    % Updates
                    tr.pi(k,l) = tr.pihat(k,l) +1/2 *p.ka(l-1)^2 *tr.w(k,l-1) *(tr.w(k,l-1) +(2 *tr.w(k,l-1) -1) *tr.da(k,l-1));

                    if tr.pi(k,l) <= 0
                        error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                    end

                    tr.mu(k,l) = tr.muhat(k,l) +1/2 *1/tr.pi(k,l) *p.ka(l-1) *tr.w(k,l-1) *tr.da(k,l-1);

                    % Volatility prediction error
                    tr.da(k,l) = (1/tr.pi(k,l) +(tr.mu(k,l) -tr.muhat(k,l))^2) *tr.pihat(k,l) -1;
                end

            elseif strcmp(r.c_prc.(type).priortype,'state') % Kalman
                % Gain update - optimal gain is calculated from ratio of input
                % variance to representation variance

                % Same gain function modified by two different alphas
                tr.g(k) = (tr.g(k-1) +tr.al(k)*p.expom)/(tr.g(k-1) +tr.al(k)*p.expom +1);
                % Hidden state mean update
                tr.mu(k,1) = tr.muhat(k,1)+tr.g(k)*tr.dau(k);
                tr.pi(k,1) = (1-tr.g(k)) *tr.al(k)*p.expom; 

                % Alternative: separate gain functions for each stimulus type
          %      if r.c_prc.(type).one_alpha
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
                tr.da(k,1) = tr.mu(k,1) -tr.muhat(k,1);

            end

            % RESPONSE BIAS
            if isfield(p,'rb')
                tr.mu(k,1) = tr.mu(k,1)+p.rb;
            end

        else
            tr.mu(k,:) = tr.mu(k-1,:); 
            tr.pi(k,:) = tr.pi(k-1,:);

            tr.muhat(k,:) = tr.muhat(k-1,:);
            tr.pihat(k,:) = tr.pihat(k-1,:);

            tr.v(k,:)  = tr.v(k-1,:);
            tr.w(k,:)  = tr.w(k-1,:);
            tr.da(k,:) = tr.da(k-1,:);
            tr.dau(k) = tr.dau(k-1);
            if isfield(tr,'g')
                tr.g(k,:)=tr.g(k-1,:);
            end
            tr.al(k)  = tr.al(k-1);
        end
        
        % Repack parameters
        for pn=1:length(pnames)
            if isfield(M.(type).p,pnames{pn})
                M.(type).p.(pnames{pn}) = p.(pnames{pn});
            end
        end

        % Repack traj
        for tn=1:length(tnames)
            if isfield(M.(type).tr,tnames{pn})
                M.(type).tr.(tnames{tn}) = tr.(tnames{tn});
            end
        end
        
        % Repack like parameters
        type='like';
        for pn=1:length(pnames)
            if isfield(M.(type).p,pnames{pn})
                M.(type).p.(pnames{pn}) = p.(pnames{pn});
            end
        end
    end
    
    % Joint prediction if more than one model
    if r.c_prc.nModels>1
        
        % Unpack like parameters
        type='like';
        pnames = fieldnames(M.(type).p);
        for pn=1:length(pnames)
            p.(pnames{pn}) = M.(type).p.(pnames{pn});
        end
        
        % set alpha
        if ~isfield(p,'al1')
            p.al1=p.al0;
        end
        if u(k,1)==0
            tr.al(k)=p.al0(u(k,2));
        elseif u(k,1)==1
            tr.al(k)=p.al1(u(k,2));
        end
        
        v_mu=[];
        v_phi=[];
        v_mu_phi=[];
        for m=1:r.c_prc.nModels

            type = r.c_prc.type{m};
            l=M.(type).l;
            
            % define mu phi
            v_mu_str = [c.(type).jointrep '(' c.(type).jointrepk ',' c.(type).jointreplev ')'];
            eval(['v_mu(m) = ' v_mu_str ';']);
            v_phi(m) = M.(type).p.phi;
            v_mu_phi(m) = v_phi(m)*v_mu(m);
        
        end
    
        % joint probability
        vj_phi = sum(v_phi);
        M.like.tr.vj_mu(k,1) = sum(v_mu_phi)/vj_phi;

        % joint prediction
        r=exp((p.eta1-M.like.tr.vj_mu(k,1))^2 - (p.eta0-M.like.tr.vj_mu(k,1))^2)/vj_phi^-2;
        M.like.tr.xchat(k,1) = 1/(1+r);
        
        % Likelihood functions: one for each
        % possible signal
        if r.c_prc.n_inputcond >1
            und1 = exp(-(u(k,1) -p.eta1)^2/(2*p.al1(u(k,2))));
            und0 = exp(-(u(k,1) -p.eta0)^2/(2*p.al0(u(k,2))));
        else
            und1 = exp(-(u(k,1) -p.eta1)^2/(2*p.al1));
            und0 = exp(-(u(k,1) -p.eta0)^2/(2*p.al0));
        end

        % Update
        M.like.tr.xc(k,1) = M.like.tr.xchat(k,1) *und1 /(xchat(k,1) *und1 +(1 -M.like.tr.xchat(k,1)) *und0);
        
        % Repack like parameters
        type='like';
        for pn=1:length(pnames)
            if isfield(M.(type).p,pnames{pn})
                M.(type).p.(pnames{pn}) = p.(pnames{pn});
            end
        end
    end
    
end

%% COMPILE RESULTS

% joint trajectories (if needed)
M.like.tr.vj_mu(1) = [];
M.like.tr.xc(1) = [];
M.like.tr.xchat(1) = [];

for m=1:r.c_prc.nModels

    % Unpack likelihood parameters
    type='like';
    pnames = fieldnames(M.(type).p);
    for pn=1:length(pnames)
        p.(pnames{pn}) = M.(type).p.(pnames{pn});
    end

    type = r.c_prc.type{m};
    l=M.(type).l;
    
    % Unpack prior parameters
    pnames = fieldnames(M.(type).p);
    for pn=1:length(pnames)
        p.(pnames{pn}) = M.(type).p.(pnames{pn});
    end

    %Unpack prior traj
    tnames = fieldnames(M.(type).tr);
    for tn=1:length(tnames)
        tr.(tnames{tn}) = M.(type).tr.(tnames{tn});
    end

    if l>1
        % Implied learning rate at the first level
        sgmmu2 = tapas_sgm(tr.mu(:,2), 1);
        lr1    = diff(sgmmu2)./tr.da(2:n,1);
        lr1(tr.da(2:n,1)==0) = 0;
    end

    % Remove representation priors
    tr.mu0(1,:)  = [];
    tr.mu(1,:)  = [];
    tr.pi(1,:)  = [];
    % Check validity of trajectories
    if any(isnan(tr.mu(:))) %|| any(isnan(pi(:)))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
    else
        % Check for implausible jumps in trajectories
        % CAB: only use first 500 trials - after that changes in precision become too small
        ntrials = min(length(tr.mu),500);
        dmu = diff(tr.mu(1:ntrials,2:end));
        dpi = diff(tr.pi(1:ntrials,2:end));
        rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
        rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

        jumpTol = 16;
        if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
            error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
            disp('Use plot for diagnosis: see within function'); % plot(abs(dpi(:,2))); hold on; plot(rmdpi(:,2),'r'); hold on; plot(jumpTol*rmdpi(:,2),'g')
        end
    end

    % Remove other dummy initial values
    tr.muhat(1,:) = [];
    tr.pihat(1,:) = [];
    tr.v(1,:)     = [];
    tr.w(1,:)     = [];
    tr.da(1,:)    = [];
    tr.dau(1)     = [];
    if isfield(tr,'g')
        tr.g(1,:)  = [];
    end
    tr.al(1)     = [];

    % Create result data structure
    traj = struct;
    
    traj.like.vj_mu = M.like.tr.vj_mu;
    traj.like.xc = M.like.tr.xc;
    traj.like.xchat = M.like.tr.xchat;

    traj.(type).mu0    = tr.mu0;
    traj.(type).mu     = tr.mu;
    traj.(type).sa     = 1./tr.pi;

    traj.(type).muhat  = tr.muhat;
    traj.(type).sahat  = 1./tr.pihat;
    traj.(type).v      = tr.v;
    traj.(type).w      = tr.w;
    traj.(type).da     = tr.da;
    traj.(type).dau    = tr.dau;
    if isfield(tr,'g')
        traj.(type).g     = tr.g;
    end

    % Updates with respect to prediction
    traj.(type).ud = tr.muhat -tr.mu;

    % Psi (precision weights on prediction errors)
    tr.psi        = NaN(n-1,l+1);
    tr.psi(:,1)   = 1./(tr.al.*tr.pi(:,1)); % dot multiply only if al is a vector
    if l>1; tr.psi(:,2)   = 1./tr.pi(:,2);end
    if l>2; tr.psi(:,3:l) = tr.pihat(:,2:l-1)./tr.pi(:,3:l);end
    traj.(type).psi   = tr.psi;

    % Epsilons (precision-weighted prediction errors)
    tr.epsi        = NaN(n-1,l);
    tr.epsi(:,1)   = tr.psi(:,1) .*tr.dau;
    if l>1; tr.epsi(:,2:l) = tr.psi(:,2:l) .*tr.da(:,1:l-1);end
    traj.(type).epsi   = tr.epsi;

    % Full learning rate (full weights on prediction errors)
    tr.wt        = NaN(n-1,l);
    if l==1; lr1=tr.psi(:,1);end
    tr.wt(:,1)   = lr1;
    if l>1; tr.wt(:,2)   = tr.psi(:,2); end
    if l>2; tr.wt(:,3:l) = 1/2 *(tr.v(:,2:l-1) *diag(p.ka(2:l-1))) .*tr.psi(:,3:l); end
    traj.(type).wt   = tr.wt;

    % Create matrices for use by the observation model
%     infStates = NaN(n-1,l,11);
%     infStates(:,:,1) = traj.muhat;
%     infStates(:,:,2) = traj.sahat;
%     infStates(:,:,3) = traj.mu;
%     infStates(:,:,4) = traj.sa;
%     infStates(:,:,5) = traj.da;
%     infStates(:,:,6) = traj.epsi;
%     infStates(:,1,7) = traj.dau;
%     infStates(:,1,8) = traj.mu0;
%     infStates(:,1,9) = traj.vj_mu;
%     infStates(:,1,10) = traj.xc;
%     infStates(:,1,11) = traj.xchat;
end

return;
