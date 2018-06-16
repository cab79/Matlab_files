function [logp, yhat, res] = logrt_softmax_binary(r, infStates, pvec)
% Calculates the log-probability of log-reaction times y (in units of log-ms) according to the
% linear log-RT model developed with Louise Marshall and Sven Bestmann
% http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002575
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2014-2016 Christoph Mathys, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

try
    l = r.c_prc.(r.c_obs.model).n_priorlevels;
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
pvec = logrt_softmax_binary_transp(r, pvec);

% Create param struct
nme=r.c_obs.pnames;
nme_gen=r.c_obs.pnames_gen;
idx=r.c_obs.priormusi;

%% SOFTMAX
if strcmp(r.c_obs.response.model, 'RT-soft') || strcmp(r.c_obs.response.model,'soft')

    % Create param struct
    type='soft';
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end
%     % Softmax parameter
%     be = exp(ptrans(8));

    % Softmax: Weed irregular trials out from inferred states, responses, and inputs
    x = r.traj.(r.c_obs.model).(r.c_obs.rep)(:,1);
    x(r.irr) = [];
    ys = r.y(:,1);
    ys(r.irr) = [];

    if size(r.u,2) == 3
        r0 = r.u(:,2);
        r0(r.irr) = [];
        r1 = r.u(:,3);
        r1(r.irr) = [];
    end

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
    logp_so(reg) = log(probc);
    yh = ys.*probc +(1-ys).*(1-probc);
    yhat(reg) = yh;
    res(reg) = (ys -yh)./sqrt(yh.*(1 -yh));
end
if strcmp(r.c_obs.response.model,'RT-soft') || strcmp(r.c_obs.response.model,'RT')

    %% RT
    % Create param struct
    type='rt';
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end
    
%     % Transform parameters to their native space
%     be0  = ptrans(1);
%     be1  = ptrans(2);
%     be2  = ptrans(3);
%     be3  = ptrans(4);
%     be4  = ptrans(5);
%     be5  = ptrans(6);
%     ze   = exp(ptrans(7));

    % logRT: Weed irregular trials out from responses and inputs
    yr = r.y(:,2); % CAB: RTs are in column 2
    yr(r.irr) = [];
    u = r.u(:,1);
    u(r.irr) = [];
    
    % check inputs RTs are logged
    if any(yr>10)
        error('input RTs must be logged and not contain any extreme outliers')
    end

    % logRT: Extract trajectories of interest from infStates
    mu1hat = r.traj.(r.c_obs.model).muhat(:,1);
    sa1hat = r.traj.(r.c_obs.model).sahat(:,1);
    da1 = r.traj.(r.c_obs.model).da(:,1);
    da2 = r.traj.(r.c_obs.model).da(:,2);
    ep2 = r.traj.(r.c_obs.model).epsi(:,2);
    ep3 = r.traj.(r.c_obs.model).epsi(:,3);
    if l>1
        mu2    = r.traj.(r.c_obs.model).(mu)(:,2);
        sa2    = r.traj.(r.c_obs.model).(sa)(:,2);
    end
    if l>2
        mu3    = r.traj.(r.c_obs.model).(mu)(:,3);
    end
    
    
    % prediction error
    % ~~~~~~~~
    da1reg = abs(da1);
    da1reg(r.irr) = [];
    ep2reg = abs(ep2);
    ep2reg(r.irr) = [];
    da2reg = da2;
    da2reg(r.irr) = [];
    ep3reg = ep3;
    ep3reg(r.irr) = [];

    % Surprise: informational
    % ~~~~~~~~
    m1hreg = mu1hat;
    m1hreg(r.irr) = [];
    poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
    surp = -log2(poo);

    % Bernoulli variance (aka irreducible uncertainty, risk) 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bernv = sa1hat;
    bernv(r.irr) = [];

    if l>1 % CAB
        % Inferential variance (aka informational or estimation uncertainty, ambiguity)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
        %sigmoid_mu2 = 1/(1+exp(-mu2)); % transform down to 1st level
        %inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
        inferv = sa2; 
        inferv(r.irr) = [];
    end

    if l>2 % CAB
        % Phasic volatility (aka environmental or unexpected uncertainty)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
        pv = exp(mu3); % transform down to 1st level
        pv(r.irr) = [];
    end

    % Calculate predicted log-reaction time
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if l>2
        logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv +be5.*da1reg +be6.*ep2reg +be7.*da2reg +be8.*ep3reg;
    elseif l>1
        logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be5.*da1reg +be6.*ep2reg;
    else
        logrt = be0 +be1.*surp +be2.*bernv +be5.*da1reg +be6.*ep2reg;
    end

    % Calculate log-probabilities for non-irregular trials
    % Note: 8*atan(1) == 2*pi (this is used to guard against
    % errors resulting from having used pi as a variable).
    logp_rt(reg) = -1/2.*log(8*atan(1).*ze) -(yr-logrt).^2./(2.*ze);
    yhat(reg) = logrt;
    res(reg) = yr-logrt;
end

if ~exist('logp_so','var')
    logp_so=0;
end
if ~exist('logp_rt','var')
    logp_rt=0;
end
logp = logp_so + logp_rt;

return;
