function [logp, yhat, res] = mismatch_softmax_binary(r, infStates, ptrans)

% CAB: Number of levels
l = r.c_prc.n_levels;

% Predictions or posteriors?
pop = 1; % Default: predictions
if r.c_obs.predorpost == 2
    pop = 3; % Alternative: posteriors
end

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
reg = ~ismember(1:n,r.irr);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

%% SOFTMAX

% Softmax parameter
be = exp(ptrans(1));

% Threshold parameter
th = ptrans(2);

% Softmax: Weed irregular trials out from inferred states, responses, and inputs
x = infStates(:,1,pop);
%x    = infStates(:,1,5); % da1
x(r.irr) = [];
mm = r.y(:,3);
mm(r.irr) = [];

if size(r.u,2) == 3
    r0 = r.u(:,2);
    r0(r.irr) = [];
    r1 = r.u(:,3);
    r1(r.irr) = [];
end

class1=1*single(mm>th);%CAB 
class2=0*single(mm<=th);%CAB 
ys=class1+class2;%CAB 

% Calculate log-probabilities for non-irregular trials
% If input matrix has only one column, assume the weight (reward value)
% of both options is equal to 1
if size(r.u,2) == 2 % CAB changed to 2 from 1, as second column now contains conditions
    % Probability of observed choice
    probc = 1./(1+exp(-be.*(2.*x-1).*(2.*ys-1)));
end
% If input matrix has three columns, the second contains the weights of
% outcome 0 and the third contains the weights of outcome 1
if size(r.u,2) == 3 % CAB changed to 3 from 2, as second column now contains conditions
    % Probability of observed choice
    probc = 1./(1+exp(-be.*(r1.*x-r0.*(1-x)).*(2.*ys-1)));
end
logp(reg) = log(probc);
yh = ys.*probc +(1-ys).*(1-probc);
yhat(reg) = yh;
res(reg) = (ys -yh)./sqrt(yh.*(1 -yh));

return;
