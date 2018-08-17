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
    
    case 1 % RT surp
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1]; % from betas 1 to 8 
    
    case 2 % 1st 2 levels of PEs (abs)
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5:9]; % from betas 1 to 8 
    
    case 3 % RT surprise + 1st 2 levels of PEs (abs)
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1 5:9]; % from betas 1 to 8 
    
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
    
    case 10 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [9]; % from betas 1 to 8 
    
    case 11 % RT 
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [10]; % from betas 1 to 8 
    
    case 12 % RT uncertainty
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1 2 3 4]; % from betas 1 to 8 
    
    case 13 % RT Prediction error
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6]; 
    
    case 14 % RT Prediction error
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [7 8]; 
    
    case 15 % RT Prediction error
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [9 10]; 
    
    case 16 % RT Prediction error
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6 7 8 9 10]; 
    
    case 17 % RT Uncertainy & PE
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [1:10]; 
    
    case 18 % RT perception
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [11]; 
    
    case 19 % PE plus Mu
    S.resp_modelspec.responses = {'RT'}; % RT (response time), Ch (choice), or both
    S.resp_modelspec.params = [5 6 7 8 9 10 11 12]; 
    
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