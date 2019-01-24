function c = GBM_config_Bayesopt(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Config structure
c = struct;

%% General settings (specific to experimental design but not to the model)
% Input intervals
% If input intervals are irregular, the last column of the input
% matrix u has to contain the interval between inputs k-1 and k
% in the k-th row, and this flag has to be set to true
c.irregular_intervals = false;
% Number of target stimuli
c.n_targets = 1;
    
%%
c.nparams =[];
c.priormus=[];
c.priorsas=[];
c.nModels = length(fieldnames(S.perc_modelspec.priormodels));
c.modelnames = fieldnames(S.perc_modelspec.priormodels);
c.n_inputcond = S.perc_modelspec.likelihood.n_inputcond;
c.response.priormodel = S.perc_modelspec.response.priormodel;
c.response.rep = S.perc_modelspec.response.rep;
c.st = [];
c.pn=0;

if ~isfield(S,'failed')
    S.failed=0;
end

%% Likelihood function (general)
switch S.perc_modelspec.likelihood.type 
    case 'binary'
        % Eta0
        c.like.eta0mu = 0;
        c.like.eta0sa = 0;

        % Eta1
        c.like.eta1mu = 1;
        c.like.eta1sa = 0;

        %% Input variance: Alpha
        switch S.perc_modelspec.likelihood.inputvar
            case 'uncertain_unequal'
                c.like.logal0mu = repmat(log(1),1,c.n_inputcond);
                c.like.logal0sa = repmat(1,1,c.n_inputcond); % unfixed
                c.like.logal0var = true; % this is a variance parameter
                c.like.logal1mu = repmat(log(1),1,c.n_inputcond);
                c.like.logal1sa = repmat(1,1,c.n_inputcond); % unfixed
                c.like.logal1var = true; % this is a variance parameter
            case 'uncertain_equal'
                % only specify al0
                c.like.logal0mu = log(0.2*ones(1,4));%repmat(log(1),1,c.n_inputcond);
                c.like.logal0sa = 100*repmat(1,1,c.n_inputcond); % unfixed
                c.like.logal0var = true; % this is a variance parameter
            case 'certain'
                % only specify al0
                c.like.logal0mu = repmat(log(0),1,c.n_inputcond); 
                c.like.logal0sa = repmat(0,1,c.n_inputcond); % fixed
                c.like.logal0var = true; % this is a variance parameter
        end
end

type='like';
c = paramvec(c,type);

%% Prior function
for m = 1:c.nModels
    
    type = c.modelnames{m};
    
    % Model name
    c.type{m} = type;

    %% add modelspecs to config struct
    for fn = fieldnames(S.perc_modelspec.priormodels.(type))'
       c.(type).(fn{1}) = S.perc_modelspec.priormodels.(type).(fn{1});
    end
    
    %% Priors
    switch S.perc_modelspec.priormodels.(type).priortype
        case 'hierarchical'
            % Initial hidden state mean and variance
                % HIERARCHICAL MODELS
                % For all but the first two levels, this is usually best
                % kept fixed to 1 (determines origin on x_i-scale). The 
                % first level is NaN because it is determined by the second,
                % and the second implies neutrality between outcomes when it
                % is centered at 0.
            
            switch S.perc_modelspec.priormodels.(type).priorupdate
                case 'dynamic'
                    switch S.perc_modelspec.priormodels.(type).n_priorlevels
                        case 1
                            c.(type).mu_0mu = [NaN, 0];
                            c.(type).mu_0sa = [NaN, 0];% prior fixed to 0 at k=0 but dynamically updates after that
                            c.(type).logsa_0mu = [NaN, log(0.5)];
                            c.(type).logsa_0sa = [NaN, 0];
                            c.(type).logsa_0var = true; % this is a variance parameter
                        case 2
                            c.(type).mu_0mu = [NaN, 0, 1];
                            c.(type).mu_0sa = [NaN, 0, 0];% prior fixed to 0 at k=0 but dynamically updates after that
                            c.(type).logsa_0mu = [NaN, log(0.5), log(1)];
                            c.(type).logsa_0sa = [NaN, 0, 0];
                            c.(type).logsa_0var = true; % this is a variance parameter
                    end

                case 'static'
                    switch S.perc_modelspec.priormodels.(type).n_priorlevels
                        case 1
                            c.(type).mu_0mu = [NaN, 0];
                            c.(type).mu_0sa = [NaN, 1];% estimated 
                            c.(type).logsa_0mu = [NaN, log(0.5)];
                            c.(type).logsa_0sa = [NaN, 0];
                            c.(type).logsa_0var = true; % this is a variance parameter
                        case 2
                            c.(type).mu_0mu = [NaN, 0, 1];
                            c.(type).mu_0sa = [NaN, 1, 0];% estimated mean
                            c.(type).logsa_0mu = [NaN, log(0.5), log(1)];
                            c.(type).logsa_0sa = [NaN, 0, 0];
                            c.(type).logsa_0var = true; % this is a variance parameter
                    end
                    
            end
            
            if ~(strcmp(S.perc_modelspec.priormodels.(type).priorupdate,'static') && S.perc_modelspec.priormodels.(type).n_priorlevels==1)
                % don't include these for static models with a single prior
            
                % Rhos
                % Format: row vector of length n_levels.
                % Undefined (therefore NaN) at the first level.
                % Fix this to zero to turn off drift.
                c.(type).rhomu = [NaN, repmat(0, 1, length(c.(type).mu_0mu)-1)];
                c.(type).rhosa = [NaN, repmat(0, 1, length(c.(type).mu_0mu)-1)];

                % Kappas
                % Format: row vector of length n_levels-1.
                % Undefined (therefore NaN) at the first level.
                % This should be fixed (preferably to 1) if the observation model
                % does not use mu_i+1 (kappa then determines the scaling of x_i+1).
                c.(type).logkamu = [NaN, repmat(log(1), 1, length(c.(type).mu_0mu)-2)];
                c.(type).logkasa = [NaN, repmat(0, 1, length(c.(type).mu_0mu)-2)];
                c.(type).logkavar = true; % this is a variance parameter

                % Format: row vector of length n_levels.
                % Undefined (therefore NaN) at the first level.
                c.(type).ommu = [NaN, repmat(-5.75-S.failed, 1, length(c.(type).mu_0mu)-2),-5.34];
                c.(type).omsa = [NaN, repmat(100, 1, length(c.(type).mu_0mu)-1)]; %abs(c.(type).ommu);
            end

        case 'state'
            % PRIOR STATE MODELS
            % e.g. Kalman Filter or Rascorla-Wagner
            % 99991   Value of the first input
            %         Usually a good choice for mu_0mu(1)
            % 99992   Variance of the first 20 inputs
            %         Usually a good choice for mu_0sa(1)
            c.(type).mu_0mu = 0.5;
            c.(type).mu_0sa = 1;
            
            % Initial Kalman gain
            c.(type).logg_0mu = 0.1;
            c.(type).logg_0sa = 0;
            
            % Format: row vector of length n_levels.
            % Undefined (therefore NaN) at the first level.
            c.(type).ommu = -6;
            c.(type).omsa = 4^2;
    end
    
    %% Response bias
    %c.(type).rbmu = 0;
    %c.(type).rbsa = 0;
    
    %% Joint models
    % Variance (Phi)
    if c.nModels>1
        c.(type).logphimu = log(2); % range between 1 and 3 seems to work
        c.(type).logphisa = 10; % unfixed
        c.(type).logphivar = true; % this is a variance parameter
    end

    %% Parameter names
    % only names added here will be estimated
    % second column: 1 if it is a variance parameter
%     c.(type).pnames = {
%         'eta0',0
%         'eta1',0
%         'logal0',1
%         'logal1',1
%         'mu_0',0
%         'logsa_0',1
%         'om',0
%         'rho',0
%         'logka',1
%         'rb',0
%         };

    % Gather prior settings and their indices in vectors
    c = paramvec(c,type);
end

% Check whether we have the right number of priors
%expectedLength = 2*c.n_levels+c.n_inputs*4+1;
%if length([c.priormus, c.priorsas]) ~= 2*expectedLength;
%    error('tapas:hgf:PriorDefNotMatchingLevels', 'Prior definition does not match number of levels.')
%end

%c.model = 'GBM';

%% Model function handle
c.prc_fun = @GBM;

% Handle to function that transforms perceptual parameters to their native space
% from the space they are estimated in
c.transp_prc_fun = @GBM_transp;

return;

function c = paramvec(c,type)
fn=fieldnames(c.(type));
for i = 1:length(fn)
    if strcmp(fn{i}(end-1:end),'mu')
        c.pn=c.pn+1;
        c.pnames{c.pn,1} = [type '_' fn{i}(1:end-2)];
        nme_gen = strsplit(fn{i}(1:end-2),'log');
        c.pnames_gen{c.pn,1} = nme_gen{end};
        c.pnames_mod{c.pn,1} = [type '_' nme_gen{end}];
        eval(['c.priormus = [c.priormus c.(type).' fn{i} '];']);
        eval(['c.nparams(c.pn) = length(c.(type).' fn{i} ');']);
        if isfield(c.(type),[fn{i}(1:end-2) 'var'])
            c.varparam(c.pn)=1;
        else
            c.varparam(c.pn)=0;
        end
    elseif strcmp(fn{i}(end-1:end),'sa')
        eval(['c.priorsas = [c.priorsas c.(type).' fn{i} '];']);
    else
        continue
    end
    if isempty(c.st)
        c.st = 0;
    else
        c.st=sum(c.nparams(1:c.pn-1));
    end
    c.priormusi{c.pn} = c.st+1:sum(c.nparams(1:c.pn));
end
