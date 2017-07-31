function y = logrt_softmax_binary_sim(r, infStates, p)
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

% CAB: Number of levels
l = r.c_prc.n_levels;

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);

y = nan(n,1);

% Initialize random number generator
rng('shuffle');

%% SOFTMAX
if strcmp(r.c_obs.model,'RT-soft') || strcmp(r.c_obs.model,'soft')
    
    % Predictions or posteriors?
    pop = 1; % Default: predictions
    if r.c_obs.predorpost == 2
        pop = 3; % Alternative: posteriors
    end

    % Softmax parameter
    be = p(8);

    % Softmax: Weed irregular trials out from inferred states, responses, and inputs
    x = infStates(:,1,pop);

    % Calculate log-probabilities for non-irregular trials
    % If input matrix has only one column, assume the weight (reward value)
    % of both options is equal to 1
    if size(x,2) == 1
        if ~any(x<0) && ~any(x>1)
            % Apply the logistic sigmoid to the inferred states
            prob = tapas_sgm(be.*(2.*x-1),1);
        else
            error('tapas:hgf:SoftMaxBinary:InfStatesIncompatible', 'infStates incompatible with tapas_softmax_binary observation model.')
        end
    end
    % If input matrix has three columns, the second contains the weights of
    % outcome 0 and the third contains the weights of outcome 1
    if size(x,2) == 2 
        % Apply the logistic sigmoid to the inferred states
        prob = tapas_sgm(be.*(x(:,1)-x(:,2)),1);
    end
    
    % Simulate
    y(:,1) = binornd(1, prob);
end
if strcmp(r.c_obs.model,'RT-soft') || strcmp(r.c_obs.model,'RT')

    %% RT

    % Transform parameters to their native space
    be0  = p(1);
    be1  = p(2);
    be2  = p(3);
    be3  = p(4);
    be4  = p(5);
    be5  = p(6);
    ze   = p(7);

    u = r.u(:,1);

    % logRT: Extract trajectories of interest from infStates
    mu1hat = infStates(:,1,1);
    sa1hat = infStates(:,1,2);
    da1 = infStates(:,1,5);
    if l>1
        mu2    = infStates(:,2,3);
        sa2    = infStates(:,2,4);
    end
    if l>2
        mu3    = infStates(:,3,3);
    end
    
    % prediction error
    % ~~~~~~~~
    da1reg = da1;
    da1reg(r.irr) = [];

    % Surprise
    % ~~~~~~~~
    m1hreg = mu1hat;
    poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
    surp = -log2(poo);

    % Bernoulli variance (aka irreducible uncertainty, risk) 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bernv = sa1hat;

    if l>1 % CAB
        % Inferential variance (aka informational or estimation uncertainty, ambiguity)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
    end

    if l>2 % CAB
        % Phasic volatility (aka environmental or unexpected uncertainty)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pv = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3); % transform down to 1st level
    end

    % Calculate predicted log-reaction time
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    try % CAB
        logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv +be5.*da1reg;
    catch
        if l<2
            logrt = be0 +be1.*surp +be2.*bernv +be5.*da1reg;
        elseif l<3
            logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be5.*da1reg;
        end
    end
    
    % Simulate
    y(:,end) = logrt+sqrt(ze)*randn(n, 1);
end


return;
