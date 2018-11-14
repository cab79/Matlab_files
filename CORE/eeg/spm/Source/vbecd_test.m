% Parameter options:
% Use time window vs. individual estimate of exact latency
% prior location and variance

dbstop if error
clear all
Dpath = 'C:\Data\CORE\eeg\ana\spm\SPMdata';
D=spm_eeg_load(fullfile(Dpath,'mspm12_fnums_flip_CORE016_4_merged_cleaned.mat'));
S.val = 1;
S.timewin = [30 50]; % time window in ms (averages over this window so should be short. Can be a single number)
S.cond = {{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'}}; % indices of trial types to run (separate dipole for each)
S.prior(1).prior_loc = [24 -34 64]; % dipole prior location in MNI. Leave empty to have a non-informative prior.
S.prior(1).prior_var = [30 30 30]; % dipole prior location variance in MNI (mm).
S.prior(1).prior_mom = []; % prior moment
S.prior(1).prior_momvar = []; % prior moment variance
S.prior(1).sym = 1; % 1=single dipole, 2 = symmetric pair
S.Niter = 10;
D = spm_eeg_inv_vbecd_cab(D,S);

% extract data
D.inv{end}.inverse.mniloc
D.inv{end}.inverse.jmni