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
% 9: level 3 prediction error
% 10: level 3 prediction error (precision-weighted)
% 11: posterior expectation

switch S.resp_model
 
    case 1 % RT Prediction error, unweighted, real
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 7 9]; 
    S.resp_modelspec.PE_abs = 0;
    
    case 2 % RT Prediction error, unweighted, abs
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 7 9]; 
    S.resp_modelspec.PE_abs = 1;
    
    case 3 % RT Prediction error, weighted, real
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [6 8 10]; 
    S.resp_modelspec.PE_abs = 0;
    
    case 4 % RT Prediction error, weighted, abs
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [6 8 10]; 
    S.resp_modelspec.PE_abs = 1;
    
    case 5 % RT uncertainty, 3 para
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [2 3 4]; 
    
    case 6 % RT uncertainty, 4 para
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1 2 3 4]; 
    
    case 7 % RT Mu
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [12 13 14]; 
    
    
    %% eeg
    
    case 22 % uncertainty
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1 2 3 4]; % from betas 1 to 8 
    
    case 23 % Prediction error
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6]; 
    
    case 24 % Prediction error
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [7 8]; 
    
    case 25 % Prediction error
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [9 10]; 
    
    case 26 % Prediction error
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6 7 8 9 10]; 
    
    case 27 % Uncertainy & PE
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1:10]; 
    
    case 28 % perception
    S.resp_modelspec.responses = {'EEG'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [11]; 
    
    case 36
    S.resp_modelspec.responses = {'EEG'};
    S.resp_modelspec.params = [5]; 
    
    case 37 
    S.resp_modelspec.responses = {'EEG'}; 
    S.resp_modelspec.params = [6]; 
    
end