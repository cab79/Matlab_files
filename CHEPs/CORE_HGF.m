cd('m:\Matlab\CORE\Pilot data\CORE_Gill2');
clear all

% load inputs
%u = load('example_cont_input.txt');
load('dt_test gi2_part2_20160727T111736.mat');
load('RT_test gi2_part2_20160727T111736.mat');
u = dt.design(3,:)';
comb = [u';RT(1,:)];

% event numbers for each condition (1st = change, 2nd = no change)
cond1 = [1 3]; %0.1 prob, 1 digit change, left hand
cond2 = [2 4]; %0.1 prob, 3 digit change, left hand
cond3 = [5 7]; %0.1 prob, 1 digit change, right hand
cond4 = [6 8]; %0.1 prob, 3 digit change, right hand
cond5 = [9 11]; %0.3 prob, 1 digit change, left hand
cond6 = [10 12]; %0.3 prob, 3 digit change, left hand
cond7 = [13 15]; %0.3 prob, 1 digit change, right hand
cond8 = [14 16]; %0.3 prob, 3 digit change, right hand
cond9 = [17 19]; %0.5 prob, 1 digit change, left hand
cond10 = [18 20]; %0.5 prob, 3 digit change, left hand
cond11 = [21 23]; %0.5 prob, 1 digit change, right hand
cond12 = [22 24]; %0.5 prob, 3 digit change, right hand

ISI = 1.0;
thresh = 0.2; % min RT that is realistic
num_trial_lag = [1 1 1]; % number of trials over which response is allowed to lag for each of 0.1 prob, 0.3 prob, 0.5 prob.

for i = 1:size(design,2)-1
    if any(cond1 == design(2,i)) || any(cond2 == design(2,i)) || any(cond3 == design(2,i)) || any(cond4 == design(2,i)) 
        lag = num_trial_lag(1);
    elseif any(cond5 == design(2,i)) || any(cond6 == design(2,i)) || any(cond7 == design(2,i)) || any(cond8 == design(2,i)) 
        lag = num_trial_lag(2);
    elseif any(cond9 == design(2,i)) || any(cond10 == design(2,i)) || any(cond11 == design(2,i)) || any(cond12 == design(2,i)) 
        lag = num_trial_lag(3);
    end
    comb(3,i)=0;
    if comb(1,i)>0
        for j = 1:lag+1
            if comb(3,i)==0 && comb(2,i+(j-1))>0 && comb(3,i-1)~=ISI*1+comb(2,i+(j-1)) && comb(2,i+(j-1))+(j-1)*ISI>thresh
                comb(3,i)=comb(2,i+(j-1)) + (j-1)*ISI;
            end
        end
    end
    if comb(1,i)==0 && comb(2,i)>0
        if ~any(comb(1,i-(j-1):i)>0)
            if comb(2,i)>thresh;
                comb(3,i)=comb(2,i);
            elseif comb(3,i-1)==0
                comb(3,i-1)=comb(2,i)+1*ISI;
            end
        end
    end
end

%conds_ana = [1 3 5];
conds_ana = [2 4 6];
idx=[];
for i = 1:length(conds_ana)
    eval(['conds = cond' num2str(conds_ana(i)) ';']); 
    for j = 1:length(conds)
        idx = [idx find(design(2,:)==conds(j))];
    end
end
idx = sort(idx);
design = design(:,idx);
comb = comb(:,idx);  
y = comb(3,:)>0; 
y=y';
u = comb(1,:)>0;
u=u';

%for i = 1:length(u)
%    if u(i)==0.1
%        u(i)=1.75;
%    end
%end

% find the Bayes optimal perceptual parameters
% prc: perceptual; obs:observation; opt:optimisation
prc_model = 'tapas_hgf_config';
obs_model = 'tapas_softmax_binary_config';
opt_algo = 'tapas_quasinewton_optim_config';
bopars = tapas_fitModel(y, u, prc_model, obs_model, opt_algo);
%tapas_fit_plotCorr(bopars)
tapas_hgf_plotTraj(bopars)

%alpha = 10; %uncertain
%alpha = 0.5; %a bit uncertain
alpha = 0.05; %certain
sim = tapas_simModel(u(2:end), 'tapas_hgf', [-0.0197 1 1 0.09 1 1 0 0 0 1 1 -1.1537 -6.2314 -6.2823 alpha], 'tapas_softmax_binary',[log(48)]);
tapas_hgf_plotTraj(sim)

%signal detection values
input = u;
%response=sim.y;
%[d,beta,C] = signal_detection(input(1:end-1),response(2:end))
response=y;
[d,beta,C] = signal_detection(input,response)

% change the responses simulated to make 0/1 changes less perceptible
prc_model = 'tapas_hgf_config';
%obs_model = 'tapas_gaussian_obs_config';
obs_model = 'tapas_softmax_binary_config';
opt_algo = 'tapas_quasinewton_optim_config';
bopars = tapas_fitModel(sim.y, u, prc_model, obs_model, opt_algo);
%bopars = tapas_fitModel(zeros(length(sim.y),1), u, prc_model, obs_model, opt_algo);
tapas_hgf_plotTraj(bopars)

% change responses to 

% bayes optimal
prc_model = 'tapas_hgf_config';
obs_model = 'tapas_bayes_optimal_binary_config';
opt_algo = 'tapas_quasinewton_optim_config';
bopars = tapas_fitModel([], u, prc_model, obs_model, opt_algo);
%tapas_fit_plotCorr(bopars)
tapas_hgf_binary_plotTraj(bopars)
