function c = logrt_softmax_binary_softmu0_config(r)
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

% Config structure
c = struct;

% Perceptual model used
c.model = r.c_prc.response.priormodel;
c.response.model = r.c_prc.response.model;

% CAB: Number of levels
try
    l = r.c_prc.(c.model).n_priorlevels;
catch
    l=1;
end

% Decision based on which representation?
c.rep = r.c_prc.response.rep; 

c.nparams =[];
c.priormus=[];
c.priorsas=[];
c.st = [];
c.pn=0;

if strcmp(c.response.model, 'RT-soft') || strcmp(c.response.model,'soft')
    % Beta
    c.soft.logbemu = log(48);
    c.soft.logbesa = 1;
    
    % Gather prior settings in vectors
    type='soft';
    c = paramvec(c,type);
end

if strcmp(c.response.model,'RT-soft') || strcmp(c.response.model,'RT')
    % Sufficient statistics of Gaussian parameter priors
    %
    % Beta_0
    c.rt.be0mu = log(500); 
    c.rt.be0sa = 4;

    % Beta_1
    c.rt.be1mu = 0;
    c.rt.be1sa = 4;

    % Beta_2
    c.rt.be2mu = 0; 
    c.rt.be2sa = 4;

    % Beta_3
    c.rt.be3mu = 0; 
    if l>1
        c.rt.be3sa = 4;
    else
        c.rt.be3sa = 0;
    end

    % Beta_4
    c.rt.be4mu = 0; 
    if l>2
        c.rt.be4sa = 4;
    else
        c.rt.be4sa = 0;
    end

    % Beta_5
    c.rt.be5mu = 0; 
    c.rt.be5sa = 4;

    % Zeta
    c.rt.logzemu = log(log(20));
    c.rt.logzesa = log(2);
    c.rt.logzesavar = true; % this is a variance parameter
    
    % Gather prior settings in vectors
    type='rt';
    c = paramvec(c,type);
end

% Model filehandle
c.obs_fun = @logrt_softmax_binary;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @logrt_softmax_binary_transp;

return;

function c = paramvec(c,type)
fn=fieldnames(c.(type));
for i = 1:length(fn)
    if strcmp(fn{i}(end-1:end),'mu')
        c.pn=c.pn+1;
        c.pnames{c.pn,1} = [type '_' fn{i}(1:end-2)];
        nme_gen = strsplit(fn{i}(1:end-2),'log');
        c.pnames_gen{c.pn,1} = nme_gen{end};
        c.pnames_mod{c.pn,1} = [type '_' nme_gen{end}];
        eval(['c.priormus = [c.priormus c.(type).' fn{i} '];']);
        eval(['c.nparams(c.pn) = length(c.(type).' fn{i} ');']);
        if isfield(c.(type),[fn{i}(1:end-2) 'var'])
            c.varparam(c.pn)=1;
        else
            c.varparam(c.pn)=0;
        end
    elseif strcmp(fn{i}(end-1:end),'sa')
        eval(['c.priorsas = [c.priorsas c.(type).' fn{i} '];']);
    else
        continue
    end
    if isempty(c.st)
        c.st = 0;
    else
        c.st=sum(c.nparams(1:c.pn-1));
    end
    c.priormusi{c.pn} = c.st+1:sum(c.nparams(1:c.pn));
end
