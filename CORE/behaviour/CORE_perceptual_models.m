function S=CORE_perceptual_models(S)

S.perc_modelspec = struct;

switch S.perc_model
    
    %% STATIC MODELS
    
%     case 1 % static prior (1 level), no sensory uncertainty
%     S.perc_modelspec.likelihood.type = 'binary'; 
%     S.perc_modelspec.likelihood.inputvar = 'certain'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 1; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'static'; % static, dynamic (unique estimate on each trial)
%     
%     case 2 % static prior (1 level), with equal sensory uncertainty for each input
%     S.perc_modelspec.likelihood.type = 'binary'; 
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 1; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'static'; % static, dynamic (unique estimate on each trial)
%     
%     case 3 % static prior (1 level), with unequal sensory uncertainty for each input
%     S.perc_modelspec.likelihood.type = 'binary'; 
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 1; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'static'; % static, dynamic (unique estimate on each trial)
    
    %% DYNAMIC MODELS: one alpha for all conds
    
%     case 1
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PR'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     case 2
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PR'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
    case 1
    S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
    S.perc_modelspec.likelihood.inputvar = 'uncertain'; % uncertain (variance), certain
    S.perc_modelspec.likelihood.inputfactors = []; % Number of conditions with unique input variance
    S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
    S.perc_modelspec.response.rep = 'mu0';
    modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
    S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
    S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
    S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
%     case 4
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     case 5
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = {'like','PL'}; % Which model representations are used for the response model? 1st: Choice. 2nd: RT
%     S.perc_modelspec.response.rep = 'xc';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     case 6
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = {'like','PL'}; % Which model representations are used for the response model? 1st: Choice. 2nd: RT
%     S.perc_modelspec.response.rep = 'xc';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     case 7
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = {'like','PR'}; % Which model representations are used for the response model? 1st: Choice. 2nd: RT
%     S.perc_modelspec.response.rep = 'xc';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     case 8
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_unequal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 1; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = {'like','PR'}; % Which model representations are used for the response model? 1st: Choice. 2nd: RT
%     S.perc_modelspec.response.rep = 'xc';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
%     
%     %% DYNAMIC MODELS: one alpha per DC
%     case 9
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 2; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PR'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
    case 2
    S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
    S.perc_modelspec.likelihood.inputvar = 'uncertain'; % uncertain (variance), certain
    S.perc_modelspec.likelihood.inputfactors = [2]; % Number of conditions with unique input variance
    S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
    S.perc_modelspec.response.rep = 'mu0';
    modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
    S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
    S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
    S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
    
    %% DYNAMIC MODELS: one alpha per DC and per side
%     case 11
%     S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
%     S.perc_modelspec.likelihood.inputvar = 'uncertain_equal'; % uncertain_equal (variance), uncertain_unequal, certain
%     S.perc_modelspec.likelihood.n_inputcond = 4; % Number of conditions with unique input variance
%     S.perc_modelspec.response.priormodel = 'PR'; % Which model representations are used for the response model?
%     S.perc_modelspec.response.rep = 'mu0';
%     modeltype='PR'; % AL (associative learning), PL (perceptual learning), PR (priming)
%     S.perc_modelspec.priormodels.(modeltype).priortype = 'state'; % constant, hierarchical, state
%     S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 0; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
%     S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
    case 3
    S.perc_modelspec.likelihood.type = 'binary'; % binary, continuous
    S.perc_modelspec.likelihood.inputvar = 'uncertain'; % uncertain (variance), certain
    S.perc_modelspec.likelihood.inputfactors = [1 2]; % Number of conditions with unique input variance
    S.perc_modelspec.response.priormodel = 'PL'; % Which model representations are used for the response model?
    S.perc_modelspec.response.rep = 'mu0';
    modeltype='PL'; % AL (associative learning), PL (perceptual learning), PR (priming)
    S.perc_modelspec.priormodels.(modeltype).priortype = 'hierarchical'; % constant, hierarchical, state
    S.perc_modelspec.priormodels.(modeltype).n_priorlevels = 2; % in prior hierarchy. For binary models, 2 is minimum; for continuous 1 is minimum.
    S.perc_modelspec.priormodels.(modeltype).priorupdate = 'dynamic'; % static, dynamic (unique estimate on each trial)
    
end