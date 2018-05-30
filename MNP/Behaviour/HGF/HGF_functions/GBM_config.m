function c = GBM_config(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Config structure
c = struct;
% Model name
c.model = 'GBM';

%% General settings (specific to experimental design but not to the model)
% Input intervals
% If input intervals are irregular, the last column of the input
% matrix u has to contain the interval between inputs k-1 and k
% in the k-th row, and this flag has to be set to true
c.irregular_intervals = false;
% Number of target stimuli
c.n_targets = 1;

    % modelspecs (add these to calling function)
    m=1;
    S.modelspec{m}.type = 'AL'; % AL (associative learning), PL (perceptual learning), PR (priming)
    S.modelspec{m}.likelihood = 'binary'; % binary, continuous
    S.modelspec{m}.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
    S.modelspec{m}.priortype = 'hierarchical'; % constant, hierarchical, state
    S.modelspec{m}.n_priorlevels = 2; % in prior hierarchy
    S.modelspec{m}.priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
%% Model settings
for m = 1:length(S.modelspec)

    % add modelspecs to config struct
    for fn = fieldnames(S.modelspec)'
       c{m}.(fn{1}) = S.modelspec{m}.(fn{1});
    end
    
    switch S.modelspec{m}.likelihood 
        case 'binary'
            % Eta0
            c.eta0mu = 0;
            c.eta0sa = 0;

            % Eta1
            c.eta1mu = 1;
            c.eta1sa = 0;

            %% Input variance: Alpha
            switch S.modelspec{m}.inputvar
                case 'uncertain_unequal'
                    c.logal0mu = log(0.05);
                    c.logal0sa = 1; % unfixed
                    c.logal1mu = log(0.05);
                    c.logal1sa = 1; % unfixed
                case 'uncertain_equal'
                    % only specify al0
                    c.logal0mu = log(0.05);
                    c.logal0sa = 1; % unfixed
                case 'certain'
                    % only specify al0
                    c.logal0mu = log(0.05); 
                    c.logal0sa = 0; % fixed
            end
    end
    
    
    %% Priors
    switch S.modelspec{m}.priortype
        case 'hierarchical'
            % Initial hidden state mean and variance
                % HIERARCHICAL MODELS
                % For all but the first two levels, this is usually best
                % kept fixed to 1 (determines origin on x_i-scale). The 
                % first level is NaN because it is determined by the second,
                % and the second implies neutrality between outcomes when it
                % is centered at 0.
            
            switch S.modelspec{m}.priorupdate
                case 'dynamic'
                    switch S.modelspec{m}.n_priorlevels
                        case 2
                            c.mu_0mu = [NaN, 0, 1];
                            c.mu_0sa = [NaN, 0, 0];% prior fixed to 0 at k=0 but dynamically updates after that
                            c.logsa_0mu = [NaN, log(0.5), log(1)];
                            c.logsa_0sa = [NaN, 0, 0];
                        case 3
                            c.mu_0mu = [NaN, 0, 1, 1];
                            c.mu_0sa = [NaN, 0, 0, 0];% prior fixed to 0 at k=0 but dynamically updates after that
                            c.logsa_0mu = [NaN, log(0.5), log(1), log(1)];
                            c.logsa_0sa = [NaN, 0, 0, 0];
                    end

                case 'static'
                    switch S.modelspec{m}.n_priorlevels
                        case 1
                            c.mu_0mu = [NaN, 0.5];
                            c.mu_0sa = [NaN, 1];% estimated mean
                            c.logsa_0mu = [NaN, log(0.1)];
                            c.logsa_0sa = [NaN, 0];
                        case 2
                            c.mu_0mu = [NaN, 0.2, 1];
                            c.mu_0sa = [NaN, 1, 0];% estimated mean
                            c.logsa_0mu = [NaN, log(0.1), log(1)];
                            c.logsa_0sa = [NaN, 0, 0];
                        case 3
                            c.mu_0mu = [NaN, 0, 1, 1];
                            c.mu_0sa = [NaN, 1, 0, 0];% estimated mean
                            c.logsa_0mu = [NaN, log(0.1), log(1), log(1)];
                            c.logsa_0sa = [NaN, 0, 0, 0];
                    end
                    
            end
            % Format: row vector of length n_levels.
            % Undefined (therefore NaN) at the first level.
            c.ommu = [NaN, repmat(-6, 1, length(c.mu_0mu)-1)];
            c.omsa = [NaN, repmat(4^2, 1, length(c.mu_0mu)-1)];
            
            % Rhos
            % Format: row vector of length n_levels.
            % Undefined (therefore NaN) at the first level.
            % Fix this to zero to turn off drift.
            c.rhomu = [NaN, repmat(0, 1, length(c.mu_0mu)-1)];
            c.rhosa = [NaN, repmat(0, 1, length(c.mu_0mu)-1)];

            % Kappas
            % Format: row vector of length n_levels-1.
            % Undefined (therefore NaN) at the first level.
            % This should be fixed (preferably to 1) if the observation model
            % does not use mu_i+1 (kappa then determines the scaling of x_i+1).
            c.logkamu = [NaN, repmat(log(1), 1, length(c.mu_0mu)-2)];
            c.logkasa = [NaN, repmat(0, 1, length(c.mu_0mu)-2)];
            
        case 'state'
            % PRIOR STATE MODELS
            % e.g. Kalman Filter or Rascorla-Wagner
            % 99991   Value of the first input
            %         Usually a good choice for mu_0mu(1)
            % 99992   Variance of the first 20 inputs
            %         Usually a good choice for mu_0sa(1)
            c.mu_0mu = 0.5;
            c.mu_0sa = 1;
            
            % Initial Kalman gain
            c.logg_0mu = 0.1;
            c.logg_0sa = 1;
            
        case 'constant'
            % must be static
            c.mu_0mu = [NaN, 0.5];
            c.mu_0sa = [NaN, 0];% constant mean
            c.logsa_0mu = [NaN, log(0.1)];
            c.logsa_0sa = [NaN, 0];
    end
    
    %% Response bias
    c.rbmu = 0;
    c.rbsa = 0;

    %% Parameter names
    % only names added here will be estimated
    % second column: 1 if it is a variance parameter
    c.pnames = {
        'eta0',0
        'eta1',0
        'logal0',1
        'logal1',1
        'mu_0',0
        'logsa_0',1
        'om',0
        'rho',0
        'logka',1
        'rb',0
        };

    %% Gather prior settings and their indices in vectors
    nparams =[];
    c.priormus=[];
    c.priorsas=[];
    st = [];
    for pn = 1:length(c.pnames)
        pname = [c.pnames{pn,1} 'mu'];
        if ~isfield(c,pname)
            continue
        end
        eval(['c.priormus = [c.priormus c.' c.pnames{pn,1} 'mu];']);
        eval(['c.priorsas = [c.priorsas c.' c.pnames{pn,1} 'sa];']);
        eval(['nparams(pn) = length(c.' c.pnames{pn,1} 'mu);']);
        if isempty(st)
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
end

return;
