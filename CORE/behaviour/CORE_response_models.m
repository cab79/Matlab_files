function S=CORE_response_models(S)

S.resp_modelspec = struct;

% betas
% 1: surprise
% 2: uncertainty (irreducible)
% 3: uncertainty (expected)
% 4: uncertainty (environmental)
% 5: level 1 prediction error
% 6: level 1 prediction error (precision-weighted)
% 7: level 2 prediction error
% 8: level 2 prediction error (precision-weighted)

switch S.resp_model
    
    case 1 % Choice
    S.resp_modelspec.responses = {'Ch'}; % RT (response time), Ch (choice), or both
    
    case 2 % RT surprise
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1]; % from betas 1 to 8 
    
    case 3 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [2]; % from betas 1 to 8 
    
    case 4 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [3]; % from betas 1 to 8 
    
    case 5 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [4]; % from betas 1 to 8 
     
    case 6 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5]; % from betas 1 to 8 
    
    case 7 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [6]; % from betas 1 to 8 
    
    case 8 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [7]; % from betas 1 to 8 
    
    case 9 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [8]; % from betas 1 to 8 
    
    case 10 % RT uncertainty
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [2 3 4]; % from betas 1 to 8 
    
    case 11 % RT Prediction error
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6 7 8]; % from betas 1 to 8 
    
    case 12 % RT full
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1:8]; % from betas 1 to 8 
end