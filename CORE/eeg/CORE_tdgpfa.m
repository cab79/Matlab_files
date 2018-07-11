% =========
% TD-GPFA DEMO 
% ========= 
% This demo shows how we can extract latent variables using TD-GPFA, and
% compare the performance of the TD-GPFA model with GPFA using both
% cross-validated log-likelihood and reconstruction errors.
% (Lakshmanan et al., Neural Computation 2015)

% Section 1 provides an example where TD-GPFA is used to extract neural
% trajectories and delays for a specified latent dimensionality.
%
% Section 2 shows how to select the optimal dimensionality using 
% cross-validation, and to compare TD-GPFA with GPFA. 


% =====
% TIPS
% =====
% For exploratory analysis using TD-GPFA, we often run only Section 1
% below, and not Section 2 (which finds the optimal latent
% dimensionality).  This can provide substantial savings in running
% time, since running Section 2 takes roughly K times as long as Section
% 1, where K is the number of cross-validation folds.  As long as we use
% a latent dimensionality that is 'large enough' in Section 1, we can
% roughly estimate the latent dimensionality by looking at the plot
% produced by plotTDGPFAlatentsVsTime.m.  The optimal latent
% dimensionality is approximately the number of top dimensions that have
% 'meaningful' temporal structure.  For visualization purposes, this
% rough dimensionality estimate is usually sufficient.
%
% The TD-GPFA model fitting requires signficiantly more computation than
% GPFA. We found that on a 2.8GhZ Intel Xeon CPU with 50GB RAM, fitting of
% TD-GPFA parameters for a dataset containing 140 trials of 45
% neurons (with 26 time bins per trial), took ~3h for 2 latent dimensions,
% and ~10h for 6 latent dimensions. However, TDGPFA can often explain the
% data with fewer latent dimensionalities than GPFA.
clear all
toolpath='C:\Data\Matlab\NeuralTraj_v0300';
addpath(genpath(toolpath));
run(fullfile(toolpath,'startup.m'));

% Data used in Simulation 1 in Lakshmanan et al., Neural Computation 2015
%set_file = 'mat_sample/sim1_single_bump';
% datFormat is set to 'seq' for our simulated data, which is continuous
% valued (see neuralTraj.m)
datFormat = 'seq';

% Sample neural spiking activity. Neural activity was recorded in the 
% primary motor cortex of a macaque monkey performing a center-out delayed
% reaching task in the Batista lab (University of Pittsburgh).  The sample
% data consist of spike trains taken from 100 ms before to 420 ms after  
% movement onset for 45 neurons and 140 reaches to the 8 targets. (For more 
% detail see Section 3.4 in Lakshmanan et al., Neural Computation 2015).
%
%dat_file = 'mat_sample/sample2_dat.mat';
% datFormat is set to 'spikes' for 0/1 spiking activity (see neuralTraj.m)
%datFormat = 'spikes';

% =========================================== 
% 1) Basic extraction of neural trajectories
% ===========================================

% Extract latent variables using TD-GPFA
method = 'tdgpfa';
dbstop if error
% Results will be saved in mat_results/runXXX/, where XXX is runIdx. Use a
% new runIdx for each dataset.
pth='C:\Data\CORE\eeg\ana\prep\cleaned\part2';
files = dir(fullfile(pth,'*.set'));
for f = 1%:length(files)
    EEG = pop_loadset('filename',files(f).name,'filepath',pth);
    datmat = EEG.data;
end

select_times=[20 400];
dp = dsearchn(EEG.times',select_times');
dsample=4;
% reformat
[nC,nD,nT] = size(datmat);
%chansample=20;
%cind = randsample(1:nC,min(chansample,nC));
try
    load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat');
catch
    chanlocs=EEG.chanlocs;
end
chans=[1:2:20,30:2:50,80:2:100];
cind = find(ismember([chanlocs(:).urchan],chans));
for t = 1:nT
    dat.seq(t).trialId = t;
    dat.seq(t).y = downsample(datmat(cind,dp(1):dp(2),t)',dsample)';
    dat.seq(t).T = size(dat.seq(t).y,2);
end
binWidth = 1000/(EEG.srate/dsample);

runIdx = 1;
parallelize=0;
%fprintf('Reading from %s \n',dat_file);
%load(dat_file);

% Choose latent dimensionality
%switch dat_file
%    case 'mat_sample/sim1_single_bump'
%        xDim = 1;
%    case 'mat_sample/sample2_dat.mat'
        xDim = 8;
%end
% NOTE: The optimal dimensionality should be found using
%       cross-validation (Section 2) below.

% Extract latent variables using TD-GPFA
ecmeMaxIters = 10;
% NOTE: To make the demo run fast, we set ecmeMaxIters to a small value.
% When actually analyzing data using TD-GPFA, we typically set ecmeMaxIters
% >= 100. Alternatively, we could omit this argument and the default value 
% of 100 will be used.
result_tdgpfa = neuralTraj(runIdx, dat, 'method', method, 'xDims', xDim,...
    'numFolds', 0, 'datFormat', datFormat, 'ecmeMaxIters',ecmeMaxIters,'parallelize',parallelize);
% NOTE: This function does most of the heavy lifting.
result_tdgpfa.binWidth = binWidth;

% Plot estimated delays
plotDelayMatrix(result_tdgpfa.estParams, result_tdgpfa.binWidth);

% Plot each dimension of neural trajectories versus time
% If sample2_dat.mat is used, xDim is set to 5 and ecmeMaxIters is set to 200, 
% for each unique target, the following plot will be similar to a panel of  
% Figure 6 in Lakshmanan et al., Neural Computation 2015.
% If sim1_single_bump.mat is used, the following plot will be similar to
% Figure 3c in Lakshmanan et al., Neural Computation 2015.
plotTDGPFAlatentsVsTime(result_tdgpfa.seqTrain, result_tdgpfa.estParams, ...
    result_tdgpfa.binWidth);

fprintf('\n');
fprintf('Basic extraction and plotting of latent variables is complete.\n');
fprintf('Press any key to start cross-validation...\n');
fprintf('[Depending on the dataset, this can take many minutes to hours.]\n');
pause;

% ========================================================
% 2) Full cross-validation to find optimal state dimensionality and compare
% TD-GPFA and GPFA models. This is analogous to Figure 5a in Lakshmanan et
% al., Neural Computation, 2015.

% Select number of cross-validation folds
numFolds = 4;

% If parallelize is true, all folds will be run in parallel using Matlab's
% parfor construct. If you have access to multiple cores, this provides
% significant speedup. NOTE: This setup works for Matlab2013b and later.
% To benefit from parallel processing, Matlab's Computing Toolbox must be
% installed. This code uses the default cluster profile.
% (see https://www.mathworks.com/help/distcomp/parallel.defaultclusterprofile.html)
parallelize = true;

% Perform cross-validation for different state dimensionalities using
% TD-GPFA. Results are saved in mat_results/runXXX/, where XXX is runIdx.
xDims = [1:8];
method = 'tdgpfa';
ecmeMaxIters = 5;
neuralTraj(runIdx, dat, 'method', method, 'xDims', xDims, 'numFolds',...
    numFolds, 'datFormat', datFormat, 'ecmeMaxIters', ecmeMaxIters, ... 
    'binWidth',binWidth,...
    'parallelize', parallelize);

% Perform cross-validation for different state dimensionalities using
% GPFA. Results are saved in mat_results/runXXX/, where XXX is runIdx.
method = 'gpfa';
emMaxIters = 5;
neuralTraj(runIdx,dat,'method',method, 'xDims', xDims,'numFolds',...
    numFolds, 'datFormat', datFormat, 'emMaxIters', emMaxIters,...
    'binWidth',binWidth,...
    'parallelize', parallelize);

% Plot cross-validated log-likelihood curves for both methods.
plotLLvsDim(runIdx);

% Compare reconstruction errors of TD-GPFA versus GPFA, analogous to Figure
% 5d in Lakshmanan et al., Neural Computation, 2015. Also plot example
% reconstructions for neurons 1,2 and 3. This is analogous to Figure 5c in
% Lakshmanan et al., Neural Computation, 2015.
neuronIds  = [1 2 3 4]; % Plot reconstructions for neurons 1, 2, 3
xDimGPFA   = 1;
xDimTDGPFA = 1;
plotReconstructions(runIdx, xDimTDGPFA, xDimGPFA, ...
    'neuronIds', neuronIds);



