function [logp, yhat, res] = mismatch_response_binary(r, infStates, ptrans)

% CAB: Number of levels
l = r.c_prc.n_levels;

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
reg = ~ismember(1:n,r.irr);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Transform parameters to their native space
be0  = ptrans(1);
be1  = ptrans(2);
be2  = ptrans(3);
be3  = ptrans(4);
be4  = ptrans(5);
be5  = ptrans(6);
ze   = exp(ptrans(7));
th   = ptrans(8);

yr = r.y(:,3); % CAB: mismatch data in column 3

% apply threshold
yr_sd = th*nanstd(yr);
yr(intersect(find(yr<yr_sd),find(yr>=-yr_sd))) = NaN;

% Weed irregular trials out from responses and inputs
yr(r.irr) = [];
u = r.u(:,1);
u(r.irr) = [];

% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
sa1hat = infStates(:,1,2);
da1    = infStates(:,1,5);
%epsi1  = infStates(:,1,6);
dau    = infStates(:,1,7);
if l>1
    mu2    = infStates(:,2,3);
    sa2    = infStates(:,2,4);
    da2    = infStates(:,2,5);
    %epsi2  = infStates(:,2,6);
end
if l>2
    mu3    = infStates(:,3,3);
    sa3    = infStates(:,3,4);
    da3    = infStates(:,3,5);
    %epsi3  = infStates(:,3,6);
end

% Surprise
% ~~~~~~~~
m1hreg = mu1hat;
m1hreg(r.irr) = [];
poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
surp = -log2(poo);

% prediction error
% ~~~~~~~~
da1reg = da1;
da1reg(r.irr) = [];
daureg = dau;
daureg(r.irr) = [];
%epsi1reg = epsi1;
%epsi1reg(r.irr) = [];
%epsi2reg = epsi2;
%epsi2reg(r.irr) = [];

if l>1 % CAB
    da2reg = da2;
    da2reg(r.irr) = [];
%    epsi3reg = epsi3;%
%	 epsi3(r.irr) = [];
end

if l>2 % CAB
    da3reg = da3;
    da3reg(r.irr) = [];
end

% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% REDUCED MODEL
logmm = be0 + be2.*daureg +be3.*da1reg;

% FULL MODEL
%try % CAB
%    logmm = be0 +be1.*surp +be2.*daureg +be3.*da1reg +be4.*da2reg +be5.*da3reg;
%catch
%    if l<2
%        logmm = be0 +be1.*surp +be2.*daureg +be3.*da1reg;
%    elseif l<3
%        logmm = be0 +be1.*surp +be2.*daureg +be3.*da1reg +be4.*da2reg;
%    end
%end

% Calculate log-probabilities for non-irregular trials
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
logp(reg) = -1/2.*log(8*atan(1).*ze) -(yr-logmm).^2./(2.*ze);
yhat(reg) = logmm;
res(reg) = yr-logmm;


return;
