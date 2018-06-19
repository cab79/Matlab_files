function [logp, yhat, res] = EEG_response_models(r, infStates, pvec)

try
    l = r.c_prc.(r.c_obs.model).n_priorlevels+1;
catch
    l=1;
end

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(r.traj.(r.c_obs.model).(r.c_obs.rep),1);
reg = ~ismember(1:n,r.irr);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Transform parameters to their native space
pvec = EEG_response_models_transp(r, pvec);

% Create param struct
nme=r.c_obs.pnames;
nme_gen=r.c_obs.pnames_gen;
idx=r.c_obs.priormusi;

type='eeg';
for pn=1:length(nme)
    if strcmp(nme{pn,1}(1:length(type)),type)
        eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
    end
end

% Weed irregular trials out from responses and inputs
ye = r.y(:,3:end); % CAB: EEG is in columns 3:end
ye(r.irr,:) = [];
u = r.u(:,1);
u(r.irr) = [];

% Extract trajectories of interest from infStates
mu1hat = r.traj.(r.c_obs.model).muhat(:,1);
sa1hat = r.traj.(r.c_obs.model).sahat(:,1);
da1 = r.traj.(r.c_obs.model).da(:,1);
da2 = r.traj.(r.c_obs.model).da(:,2);
ep2 = r.traj.(r.c_obs.model).epsi(:,2);
if l>1
    mu2    = r.traj.(r.c_obs.model).mu(:,2);
    sa2    = r.traj.(r.c_obs.model).sa(:,2);
end
if l>2
    mu3    = r.traj.(r.c_obs.model).mu(:,3);
    ep3 = r.traj.(r.c_obs.model).epsi(:,3);
end


% prediction error
% ~~~~~~~~
da1reg = abs(da1);
da1reg(r.irr) = [];
ep2reg = abs(ep2);
ep2reg(r.irr) = [];
da2reg = da2;
da2reg(r.irr) = [];
if l>2
    ep3reg = ep3;
    ep3reg(r.irr) = [];
end

% Bernoulli variance (aka irreducible uncertainty, risk) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bernv = sa1hat;
bernv(r.irr) = [];

% Precision-weighted prediction of level 1
% ~~~~~~~~
m1hreg = mu1hat;
m1hreg(r.irr) = [];
pwp1 = m1hreg.*(1./bernv); 

% level 2 prediction
m2reg = mu2;
m2reg(r.irr) = [];

% level 3 prediction
m3reg = mu3;
m3reg(r.irr) = [];

% if l>1 % CAB
%     % Inferential variance (aka informational or estimation uncertainty, ambiguity)
%     inferv = sa2; 
%     inferv(r.irr) = [];
% end
% 
% if l>2 % CAB
%     % Phasic volatility (aka environmental or unexpected uncertainty)
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     %pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
%     pv = exp(mu3); 
%     pv(r.irr) = [];
% end

% Calculate predicted log-EEG
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if l>2
    logeeg = be0 +be1.*m1hreg +be2.*pwp1 +be3.*m2reg +be4.*m3reg +be5.*da1reg +be6.*ep2reg +be7.*da2reg +be8.*ep3reg;
elseif l>1
    logeeg = be0 +be1.*m1hreg +be2.*pwp1 +be3.*m2reg +be5.*da1reg +be6.*ep2reg;
else
    logeeg = be0 +be1.*m1hreg +be2.*pwp1 +be5.*da1reg +be6.*ep2reg;
end

% Calculate log-probabilities for non-irregular trials
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
logp(reg) = -1/2.*log(8*atan(1).*ze) -(ye-logeeg).^2./(2.*ze);
yhat(reg) = logeeg;
res(reg) = ye-logeeg;

return;
