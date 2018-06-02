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
        %if r.c_prc.(type).one_alpha
        %    g  = NaN(n,1); % Kalman gain (optional)
        %else
            M.(type).tr.g  = NaN(n,2); % Kalman gain (optional)
        %end
        M.(type).tr.g(1,:)  = M.(type).p.g_0; % Kalman gain (optional)
        M.(type).p.expom = exp(M.(type).p.om);
    end
end

%% UPDATES
for k=2:1:n
    
    for m=1:r.c_prc.nModels
        if not(ismember(k-1, r.ign))
            
%             % Unpack likelihood parameters
%             type='like';
%             pnames = fieldnames(M.(type).p);
%             for pn=1:length(pnames)
%                 p.(pnames{pn}) = M.(type).p.(pnames{pn});
%             end
% 
            type = r.c_prc.type{m};
            l=M.(type).l;
%             
%             % Unpack prior parameters
%             pnames = fieldnames(M.(type).p);
%             for pn=1:length(pnames)
%                 p.(pnames{pn}) = M.(type).p.(pnames{pn});
%             end
%             
%             %Unpack prior traj
%             tnames = fieldnames(M.(type).tr);
%             for tn=1:length(tnames)
%                 tr.(tnames{tn}) = M.(type).tr.(tnames{tn});
%             end
        
            %% Predictions (from previous trial or fixed parameters)
            if strcmp(r.c_prc.(type).priortype,'hierarchical')
                if strcmp(r.c_prc.(type).priorupdate,'fixed')
                    if r.c_prc.(type).n_muin>1
                        M.(type).tr.muhat(1,2) = M.(type).p.mu_0(1+u(k,2));
                        M.(type).tr.pihat(1,2) = 1./M.(type).p.sa_0(1+u(k,2));
                    else
                        M.(type).tr.muhat(1,2) = M.(type).p.mu_0(2);
                        M.(type).tr.pihat(1,2) = 1./M.(type).p.sa_0(2);
                    end
                    % 2nd level prediction
                    M.(type).tr.muhat(k,2) = M.(type).tr.muhat(1,2);%mu(k-1,2) +t(k) *rho(2); % fixed to initial value - not updated on each trial

                elseif strcmp(r.c_prc.(type).priorupdate,'dynamic')
                    % 2nd level prediction
                    M.(type).tr.muhat(k,2) = M.(type).tr.mu(k-1,2) +t(k) *M.(type).p.rho(2);

                end
                % Prediction from level 2 (which can be either fixed or dynamic)
                M.(type).tr.muhat(k,1) = tapas_sgm(M.(type).tr.muhat(k,2), 1);

            elseif strcmp(r.c_prc.(type).priortype,'state')
                % Prediction from prior state, e.g. Kalman filter
                M.(type).tr.muhat(k,1) =  M.(type).tr.mu(k-1,1);

            end

            % Precision of prediction
            M.(type).tr.pihat(k,1) = 1/(M.(type).tr.muhat(k,1)*(1 -M.(type).tr.muhat(k,1)));


            %% Updates

            % Value prediction error, e.g. for Kalman filter, also known as the
            % "innovation"
            M.(type).tr.dau(k) = u(k,1) -M.(type).tr.muhat(k,1);

            % set alpha
            if ~isfield(M.like.p,'al1')
                M.like.p.al1=M.like.p.al0;
            end
            if u(k,1)==0
                M.like.tr.al(k)=M.like.p.al0(u(k,2));
            elseif u(k,1)==1
                M.like.tr.al(k)=M.like.p.al1(u(k,2));
            end

            % 
            if strcmp(r.c_prc.(type).priortype,'hierarchical')
                % Likelihood functions: one for each
                % possible signal
                if r.c_prc.n_inputcond >1
                    und1 = exp(-(u(k,1) -M.like.p.eta1)^2/(2*M.like.p.al1(u(k,2))));
                    und0 = exp(-(u(k,1) -M.like.p.eta0)^2/(2*M.like.p.al0(u(k,2))));
                else
                    und1 = exp(-(u(k,1) -M.like.p.eta1)^2/(2*M.like.p.al1));
                    und0 = exp(-(u(k,1) -M.like.p.eta0)^2/(2*M.like.p.al0));
                end


                
                if strcmp(type,'AL')
                    if u(k,3)==2
                        M.(type).tr.mu0(k,1) = M.(type).tr.muhat(k,1) *und1 /(M.(type).tr.muhat(k,1) *und1 +(1 -M.(type).tr.muhat(k,1)) *und0);
                        M.(type).tr.mu(k,1) = M.(type).tr.mu0(k,1);
                    elseif u(k,3)==1
                        M.(type).tr.mu0(k,1) = (1-M.(type).tr.muhat(k,1)) *und1 /(M.(type).tr.muhat(k,1) *und0 +(1 -M.(type).tr.muhat(k,1)) *und1);
                        M.(type).tr.mu(k,1) = 1-M.(type).tr.mu0(k,1);

                        % calculate prediction error for mu0 - muhat

                    end
                elseif strcmp(type,'PL')
                    M.(type).tr.mu(k,1) = M.(type).tr.muhat(k,1) *und1 /(M.(type).tr.muhat(k,1) *und1 +(1 -M.(type).tr.muhat(k,1)) *und0);
                    M.(type).tr.mu0(k,1) = M.(type).tr.mu(k,1);
                end


                %%
                % Representation prediction error
                M.(type).tr.da(k,1) = M.(type).tr.mu(k,1) -M.(type).tr.muhat(k,1);

                % second level predictions and precisions
                if strcmp(r.c_prc.(type).priorupdate,'fixed')
                    M.(type).tr.mu(k,2) = M.(type).tr.muhat(k,2); % for a model with higher level predictions, which are fixed
                    % At second level, assume Inf precision for a model with invariable predictions
                    M.(type).tr.pi(k,2) = Inf;
                    M.(type).tr.pihat(k,2) = Inf;

                elseif strcmp(r.c_prc.(type).priorupdate,'dynamic')
                    % Precision of prediction
                    M.(type).tr.pihat(k,2) = 1/(1/M.(type).tr.pi(k-1,2) +exp(M.(type).p.ka(2) *M.(type).tr.mu(k-1,3) +M.(type).p.om(2)));

                    % Updates
                    M.(type).tr.pi(k,2) = M.(type).tr.pihat(k,2) +1/M.(type).tr.pihat(k,1);
                    M.(type).tr.mu(k,2) = M.(type).tr.muhat(k,2) +1/M.(type).tr.pi(k,2) *M.(type).tr.da(k,1);

                    % Volatility prediction error
                    M.(type).tr.da(k,2) = (1/M.(type).tr.pi(k,2) +(M.(type).tr.mu(k,2) -M.(type).tr.muhat(k,2))^2) *M.(type).tr.pihat(k,2) -1;
                end

                % Implied posterior precision at first level
                sgmmu2 = tapas_sgm(M.(type).tr.mu(k,2), 1);
                M.(type).tr.pi(k,1) = M.(type).tr.pi(k,2)/(sgmmu2*(1-sgmmu2));

                if l > 3
                    % Pass through higher levels
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for j = 3:l-1
                        % Prediction
                        M.(type).tr.muhat(k,j) = M.(type).tr.mu(k-1,j) +t(k) *M.(type).p.rho(j);

                        % Precision of prediction
                        M.(type).tr.pihat(k,j) = 1/(1/M.(type).tr.pi(k-1,j) +t(k) *exp(M.(type).p.ka(j) *M.(type).tr.mu(k-1,j+1) +M.(type).p.om(j)));

                        % Weighting factor
                        M.(type).tr.v(k,j-1) = t(k) *exp(M.(type).p.ka(j-1) *M.(type).tr.mu(k-1,j) +M.(type).p.om(j-1));
                        M.(type).tr.w(k,j-1) = M.(type).tr.v(k,j-1) *M.(type).tr.pihat(k,j-1);

                        % Updates
                        M.(type).tr.pi(k,j) = M.(type).tr.pihat(k,j) +1/2 *M.(type).p.ka(j-1)^2 *M.(type).tr.w(k,j-1) *(M.(type).tr.w(k,j-1) +(2 *M.(type).tr.w(k,j-1) -1) *M.(type).tr.da(k,j-1));

                        if M.(type).tr.pi(k,j) <= 0
                            error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                        end

                        M.(type).tr.mu(k,j) = M.(type).tr.muhat(k,j) +1/2 *1/M.(type).tr.pi(k,j) *M.(type).p.ka(j-1) *M.(type).tr.w(k,j-1) *M.(type).tr.da(k,j-1);

                        % Volatility prediction error
                        M.(type).tr.da(k,j) = (1/M.(type).tr.pi(k,j) +(M.(type).tr.mu(k,j) -M.(type).tr.muhat(k,j))^2) *M.(type).tr.pihat(k,j) -1;
                    end
                end
                if l>2
                    % Last level
                    % ~~~~~~~~~~
                    % Prediction
                    M.(type).tr.muhat(k,l) = M.(type).tr.mu(k-1,l) +t(k) *M.(type).p.rho(l);

                    % Precision of prediction
                    M.(type).tr.pihat(k,l) = 1/(1/M.(type).tr.pi(k-1,l) +t(k) *M.(type).p.th);

                    % Weighting factor
                    M.(type).tr.v(k,l)   = t(k) *M.(type).p.th;
                    M.(type).tr.v(k,l-1) = t(k) *exp(M.(type).p.ka(l-1) *M.(type).tr.mu(k-1,l) +M.(type).p.om(l-1));
                    M.(type).tr.w(k,l-1) = M.(type).tr.v(k,l-1) *M.(type).tr.pihat(k,l-1);

                    % Updates
                    M.(type).tr.pi(k,l) = M.(type).tr.pihat(k,l) +1/2 *M.(type).p.ka(l-1)^2 *M.(type).tr.w(k,l-1) *(M.(type).tr.w(k,l-1) +(2 *M.(type).tr.w(k,l-1) -1) *M.(type).tr.da(k,l-1));

                    if M.(type).tr.pi(k,l) <= 0
                        error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                    end

                    M.(type).tr.mu(k,l) = M.(type).tr.muhat(k,l) +1/2 *1/M.(type).tr.pi(k,l) *M.(type).p.ka(l-1) *M.(type).tr.w(k,l-1) *M.(type).tr.da(k,l-1);

                    % Volatility prediction error
                    M.(type).tr.da(k,l) = (1/M.(type).tr.pi(k,l) +(M.(type).tr.mu(k,l) -M.(type).tr.muhat(k,l))^2) *M.(type).tr.pihat(k,l) -1;
                end

            elseif strcmp(r.c_prc.(type).priortype,'state') % Kalman
                % Gain update - optimal gain is calculated from ratio of input
                % variance to representation variance

                % Same gain function modified by two different alphas
                M.(type).tr.g(k) = (M.(type).tr.g(k-1) +M.(type).tr.al(k)*M.(type).p.expom)/(M.(type).tr.g(k-1) +M.(type).tr.al(k)*M.(type).p.expom +1);
                % Hidden state mean update
                M.(type).tr.mu(k,1) = M.(type).tr.muhat(k,1)+M.(type).tr.g(k)*M.(type).tr.dau(k);
                M.(type).tr.mu0(k,1) = M.(type).tr.mu(k,1);
                M.(type).tr.pi(k,1) = (1-M.(type).tr.g(k)) *M.like.tr.al(k)*M.(type).p.expom; 

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
                M.(type).tr.da(k,1) = M.(type).tr.mu(k,1) -M.(type).tr.muhat(k,1);

            end

            % RESPONSE BIAS
            if isfield(M.(type).p,'rb')
                M.(type).tr.mu(k,1) = M.(type).tr.mu(k,1)+M.(type).p.rb;
            end

        else
            M.(type).tr.mu(k,:) = M.(type).tr.mu(k-1,:); 
            M.(type).tr.pi(k,:) = M.(type).tr.pi(k-1,:);

            M.(type).tr.muhat(k,:) = M.(type).tr.muhat(k-1,:);
            M.(type).tr.pihat(k,:) = M.(type).tr.pihat(k-1,:);

            M.(type).tr.v(k,:)  = M.(type).tr.v(k-1,:);
            M.(type).tr.w(k,:)  = M.(type).tr.w(k-1,:);
            M.(type).tr.da(k,:) = M.(type).tr.da(k-1,:);
            M.(type).tr.dau(k) = M.(type).tr.dau(k-1);
            if isfield(M.(type).tr,'g')
                M.(type).tr.g(k,:)=M.(type).tr.g(k-1,:);
            end
            M.like.tr.al(k)  = M.like.tr.al(k-1);
        end
        
%         % Repack parameters
%         for pn=1:length(pnames)
%             if isfield(M.(type).p,pnames{pn})
%                 M.(type).p.(pnames{pn}) = p.(pnames{pn});
%             end
%         end
% 
%         % Repack traj
%         for tn=1:length(tnames)
%             if isfield(M.(type).tr,tnames{tn})
%                 M.(type).tr.(tnames{tn}) = tr.(tnames{tn});
%             end
%         end
%         
%         % Repack like parameters
%         type='like';
%         for pn=1:length(pnames)
%             if isfield(M.(type).p,pnames{pn})
%                 M.(type).p.(pnames{pn}) = p.(pnames{pn});
%             end
%         end
    end
    
    % Joint prediction if more than one model
    if r.c_prc.nModels>1
        
%         % Unpack like parameters
%         type='like';
%         pnames = fieldnames(M.(type).p);
%         for pn=1:length(pnames)
%             p.(pnames{pn}) = M.(type).p.(pnames{pn});
%         end
        
        % set alpha
        if ~isfield(M.like.p,'al1')
            M.like.p.al1=M.like.p.al0;
        end
        if u(k,1)==0
            M.like.tr.al(k)=M.like.p.al0(u(k,2));
        elseif u(k,1)==1
            M.like.tr.al(k)=M.like.p.al1(u(k,2));
        end
        
        v_mu=[];
        v_phi=[];
        v_mu_phi=[];
        for m=1:r.c_prc.nModels

            type = r.c_prc.type{m};
            l=M.(type).l;
            
            % define mu phi
            v_mu_str = ['M.(type).tr.' r.c_prc.(type).jointrep '(' r.c_prc.(type).jointrepk ',' num2str(r.c_prc.(type).jointreplev) ')'];
            eval(['v_mu(m) = ' v_mu_str ';']);
            v_phi(m) = M.(type).p.phi;
            v_mu_phi(m) = v_phi(m)*v_mu(m);
        
        end
    
        % joint probability
        vj_phi = sum(v_phi);
        M.like.tr.vj_mu(k,1) = sum(v_mu_phi)/vj_phi;

        % joint prediction
        rt=exp((M.like.p.eta1-M.like.tr.vj_mu(k,1))^2 - (M.like.p.eta0-M.like.tr.vj_mu(k,1))^2)/vj_phi^-2;
        M.like.tr.xchat(k,1) = 1/(1+rt);
        
        % Likelihood functions: one for each
        % possible signal
        if r.c_prc.n_inputcond >1
            und1 = exp(-(u(k,1) -M.like.p.eta1)^2/(2*M.like.p.al1(u(k,2))));
            und0 = exp(-(u(k,1) -M.like.p.eta0)^2/(2*M.like.p.al0(u(k,2))));
        else
            und1 = exp(-(u(k,1) -M.like.p.eta1)^2/(2*M.like.p.al1));
            und0 = exp(-(u(k,1) -M.like.p.eta0)^2/(2*M.like.p.al0));
        end

        % Update
        M.like.tr.xc(k,1) = M.like.tr.xchat(k,1) *und1 /(M.like.tr.xchat(k,1) *und1 +(1 -M.like.tr.xchat(k,1)) *und0);
        
%         % Repack like parameters
%         type='like';
%         for pn=1:length(pnames)
%             if isfield(M.(type).p,pnames{pn})
%                 M.(type).p.(pnames{pn}) = p.(pnames{pn});
%             end
%         end
    end
    
end

%% COMPILE RESULTS

% joint trajectories (if needed)
M.like.tr.vj_mu(1) = [];
M.like.tr.xc(1) = [];
M.like.tr.xchat(1) = [];

for m=1:r.c_prc.nModels

%     % Unpack likelihood parameters
%     type='like';
%     pnames = fieldnames(M.(type).p);
%     for pn=1:length(pnames)
%         p.(pnames{pn}) = M.(type).p.(pnames{pn});
%     end

    type = r.c_prc.type{m};
    l=M.(type).l;
    
%     % Unpack prior parameters
%     pnames = fieldnames(M.(type).p);
%     for pn=1:length(pnames)
%         p.(pnames{pn}) = M.(type).p.(pnames{pn});
%     end
% 
%     %Unpack prior traj
%     tnames = fieldnames(M.(type).tr);
%     for tn=1:length(tnames)
%         tr.(tnames{tn}) = M.(type).tr.(tnames{tn});
%     end

    if l>1
        % Implied learning rate at the first level
        sgmmu2 = tapas_sgm(M.(type).tr.mu(:,2), 1);
        lr1    = diff(sgmmu2)./M.(type).tr.da(2:n,1);
        lr1(M.(type).tr.da(2:n,1)==0) = 0;
    end

    % Remove representation priors
    M.(type).tr.mu0(1,:)  = [];
    M.(type).tr.mu(1,:)  = [];
    M.(type).tr.pi(1,:)  = [];
    % Check validity of trajectories
    if any(isnan(M.(type).tr.mu(:))) %|| any(isnan(pi(:)))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
    else
        % Check for implausible jumps in trajectories
        % CAB: only use first 500 trials - after that changes in precision become too small
        ntrials = min(length(M.(type).tr.mu),500);
        dmu = diff(M.(type).tr.mu(1:ntrials,2:end));
        dpi = diff(M.(type).tr.pi(1:ntrials,2:end));
        rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
        rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

        jumpTol = 16;
        if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
            error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
            disp('Use plot for diagnosis: see within function'); % plot(abs(dpi(:,2))); hold on; plot(rmdpi(:,2),'r'); hold on; plot(jumpTol*rmdpi(:,2),'g')
        end
    end

    % Remove other dummy initial values
    M.(type).tr.muhat(1,:) = [];
    M.(type).tr.pihat(1,:) = [];
    M.(type).tr.v(1,:)     = [];
    M.(type).tr.w(1,:)     = [];
    M.(type).tr.da(1,:)    = [];
    M.(type).tr.dau(1)     = [];
    if isfield(M.(type).tr,'g')
        M.(type).tr.g(1,:)  = [];
    end
    M.like.tr.al(1)     = [];

    % Create result data structure
    traj = struct;
    
    traj.like.vj_mu = M.like.tr.vj_mu;
    traj.like.xc = M.like.tr.xc;
    traj.like.xchat = M.like.tr.xchat;

    traj.(type).mu0    = M.(type).tr.mu0;
    traj.(type).mu     = M.(type).tr.mu;
    traj.(type).sa     = 1./M.(type).tr.pi;

    traj.(type).muhat  = M.(type).tr.muhat;
    traj.(type).sahat  = 1./M.(type).tr.pihat;
    traj.(type).v      = M.(type).tr.v;
    traj.(type).w      = M.(type).tr.w;
    traj.(type).da     = M.(type).tr.da;
    traj.(type).dau    = M.(type).tr.dau;
    if isfield(M.(type).tr,'g')
        traj.(type).g     = M.(type).tr.g;
    end

    % Updates with respect to prediction
    traj.(type).ud = M.(type).tr.muhat -M.(type).tr.mu;

    % Psi (precision weights on prediction errors)
    M.(type).tr.psi        = NaN(n-1,l+1);
    M.(type).tr.psi(:,1)   = 1./(M.like.tr.al.*M.(type).tr.pi(:,1)); % dot multiply only if al is a vector
    if l>1; M.(type).tr.psi(:,2)   = 1./M.(type).tr.pi(:,2);end
    if l>2; M.(type).tr.psi(:,3:l) = M.(type).tr.pihat(:,2:l-1)./M.(type).tr.pi(:,3:l);end
    traj.(type).psi   = M.(type).tr.psi;

    % Epsilons (precision-weighted prediction errors)
    M.(type).tr.epsi        = NaN(n-1,l);
    M.(type).tr.epsi(:,1)   = M.(type).tr.psi(:,1) .*M.(type).tr.dau;
    if l>1; M.(type).tr.epsi(:,2:l) = M.(type).tr.psi(:,2:l) .*M.(type).tr.da(:,1:l-1);end
    traj.(type).epsi   = M.(type).tr.epsi;

    % Full learning rate (full weights on prediction errors)
    M.(type).tr.wt        = NaN(n-1,l);
    if l==1; lr1=M.(type).tr.psi(:,1);end
    M.(type).tr.wt(:,1)   = lr1;
    if l>1; M.(type).tr.wt(:,2)   = M.(type).tr.psi(:,2); end
    if l>2; M.(type).tr.wt(:,3:l) = 1/2 *(M.(type).tr.v(:,2:l-1) *diag(M.(type).p.ka(2:l-1))) .*M.(type).tr.psi(:,3:l); end
    traj.(type).wt   = M.(type).tr.wt;

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
