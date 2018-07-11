function c = mismatch_softmax_binary_config(r)
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

% Is the decision based on predictions or posteriors? Comment as appropriate.
%c.predorpost = 1; % Predictions
c.predorpost = 2; % Posteriors

% Model name
c.model = '';

% Beta
c.logbemu = log(48);
c.logbesa = 1;

% Threshold
c.thmu = 0;
c.thsa = 2;

% Gather prior settings in vectors
c.priormus = [
    c.logbemu,...
    c.thmu,...
         ];

c.priorsas = [
    c.logbesa,...
    c.thsa,...
         ];

% Model filehandle
c.obs_fun = @mismatch_softmax_binary;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @mismatch_softmax_binary_transp;

return;
