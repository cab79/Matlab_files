function c = logrt_softmax_binary_RTsoft_config(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the linear log-reaction time response model according to as
% developed with Louise Marshall and Sven Bestmann
% http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002575
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Gaussian noise observation model assumes that responses have a Gaussian distribution around
% the inferred mean of the relevant state. The only parameter of the model is the noise variance
% (NOT standard deviation) zeta.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2016 Christoph Mathys, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% CAB: Number of levels
l = r.c_prc.n_levels;

% Config structure
c = struct;

% Model name
%c.model = 'soft';
%c.model = 'RT';
c.model = 'RT-soft';

% Is the decision based on predictions or posteriors? Comment as appropriate.
%c.predorpost = 1; % Predictions
c.predorpost = 2; % Posteriors

% Beta
c.logbemu = log(48);
c.logbesa = 1;

% Sufficient statistics of Gaussian parameter priors
%
% Beta_0
c.be0mu = log(500); 
c.be0sa = 4;

% Beta_1
c.be1mu = 0;
c.be1sa = 4;

% Beta_2
c.be2mu = 0; 
c.be2sa = 4;

% Beta_3
c.be3mu = 0; 
if l>1
    c.be3sa = 4;
else
    c.be3sa = 0;
end

% Beta_4
c.be4mu = 0; 
if l>2
    c.be4sa = 4;
else
    c.be4sa = 0;
end

% Beta_5
c.be5mu = 0; 
c.be5sa = 4;

% Zeta
c.logzemu = log(log(20));
c.logzesa = log(2);

% Gather prior settings in vectors
c.priormus = [
    c.be0mu,...
    c.be1mu,...
    c.be2mu,...
    c.be3mu,...
    c.be4mu,...
    c.be5mu,...
    c.logzemu,...
    c.logbemu,...
         ];

c.priorsas = [
    c.be0sa,...
    c.be1sa,...
    c.be2sa,...
    c.be3sa,...
    c.be4sa,...
    c.be5sa,...
    c.logzesa,...
    c.logbesa,...
         ];

% Model filehandle
c.obs_fun = @logrt_softmax_binary;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @logrt_softmax_binary_transp;

return;
