function c = EEG_response_models_config(r)

% Config structure
c = struct;

% Perceptual model used
c.model = r.c_prc.response.priormodel;
c.response.model = r.c_prc.response.model;

% CAB: Number of levels
try
    l = r.c_prc.(c.model).n_priorlevels+1;
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

% Sufficient statistics of Gaussian parameter priors
%
% Beta_0
c.eeg.be0mu = log(0); 
c.eeg.be0sa = 4;

% Beta_1
c.eeg.be1mu = 0;
c.eeg.be1sa = 4;

% Beta_2
c.eeg.be2mu = 0; 
c.eeg.be2sa = 4;

% Beta_3
c.eeg.be3mu = 0; 
if l>1
    c.eeg.be3sa = 4;
else
    c.eeg.be3sa = 0;
end

% Beta_4
c.eeg.be4mu = 0; 
c.eeg.be7mu = 0; 
c.eeg.be8mu = 0; 
if l>2
    c.eeg.be4sa = 4;
    c.eeg.be7sa = 4;
    c.eeg.be8sa = 4;
else
    c.eeg.be4sa = 0;
    c.eeg.be7sa = 0;
    c.eeg.be8sa = 0;
end

% Beta_5
c.eeg.be5mu = 0; 
c.eeg.be5sa = 4;

% Beta_6
c.eeg.be6mu = 0; 
c.eeg.be6sa = 4;

% Zeta
c.eeg.logzemu = log(log(20));
c.eeg.logzesa = log(2);
c.eeg.logzevar = true; % this is a variance parameter

% Gather prior settings in vectors
type='eeg';
c = paramvec(c,type);

% Model filehandle
c.obs_fun = @EEG_response_models;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @EEG_response_models_transp;

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
