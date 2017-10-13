function [Ttmp, ptmp, fweptmp, fdrptmp,COPE] = perform_glm_with_snpm(edges, designMatrix, contrasts, standardise,S)
%PERFORM_GLM_WITH_SNPM - CAB modification of ROI_NETS function for FSL:
%
% [T, P, CORRP] = PERFORM_GLM_WITH_RANDOMISE(DATA, X, CONTRASTS, STANDARDISE)
%   runs a GLM of DATA = X Beta + noise. Returns T-statistics, uncorrected
%   and permutation-corrected p-values associated with the regression. 
%
% [T, P, CORRP, COPE] = ... also returns the contrast of parameter
%   estimates from the regression. 
%
% Randomise actually returns 1-p values in p and corrp.
%
% Note: randomise automatically demeans, so this is not helpful for finding
% mean effect over all subjects / sessions

nSessions = cols(edges);

[checkMe, nEVs] = size(designMatrix);
assert(checkMe == nSessions, ...
       [mfilename ':BadDesign'], ...
       'Design matrix must have as many rows as the data have columns. \n');
   
[nContrasts, checkMe] = size(contrasts);
assert(checkMe == nEVs, ...
       [mfilename ':BadContrasts'], ...
       'Contrasts must have as many columns as EVs in the design matrix. \n');

if nargin < 4 || ~exist('standardise', 'var'),
    standardise = false;
else
    assert(islogical(standardise), ...
           [mfilename ':BadStandardise'], ...
           'Standardise input must be a logical value. \n');
end%if

% check for NaNs - treat as subjects to be ignored. This is a good
% representation for missing data. 
%
% What will we do with NaNs? Maybe impute the value as the group mean?
% Maybe EM imputation?
% Instead, let's remove them outright. 
badSubjects = any(isnan(edges),1);
cleanEdges  = edges;
cleanEdges(:,badSubjects) = [];

%% Save out edges into nifti
S.data_path = S.clus_path{S.cldir};
for iS = cols(cleanEdges):-1:1,
    inputNifti{iS,1} = fullfile(S.data_path, ['S' num2str(iS) '_allavg.nii']);
    nii = make_nii(cleanEdges(:,iS)', [1 1 1], [0 0 0]); 
    save_nii(nii, inputNifti{iS,1});
end
Ci = onCleanup(@() delete(inputNifti));

%% Construct design matrix
designMatrix(badSubjects,:) = [];
if standardise, 
    % demean and variance normalise
    X = bsxfun(@rdivide, bsxfun(@minus, designMatrix, mean(designMatrix)), ...
                         std(designMatrix));
else
    X = designMatrix;
end%if
for x = 1:length(X)
    grpind(x,1) = find(X(x,:));
end

% save out
%designFile = fullfile(S.datapath, 'univariate_edge_test_design.mat');
%save_vest(X, designFile);
%Cd = onCleanup(@() delete(designFile));

%% Construct contrasts
%contrastFile = fullfile(S.datapath, 'univariate_edge_test_design.con');
%save_vest(contrasts, contrastFile);
%Cc = onCleanup(@() delete(contrastFile));

%% Run randomise
%outputNifti = fullfile(S.datapath, 'univariate_edge_test');

% call to randomise
%command = sprintf('randomise -i %s -o %s -d %s -t %s -x --norcmask', ...
%                  inputNifti, outputNifti, designFile, contrastFile);
              
%% Produce nice COPEs for each edge 
% a cope is the difference in mean Z-converted correlations between each
% group
pinvxtx = pinv(designMatrix' * designMatrix);
pinvx   = pinvxtx * designMatrix';

for iEdge = rows(edges):-1:1,
    COPE(iEdge,:) = contrasts * pinvx * edges(iEdge, :).';
end%for

%% submit to design_batch
% specify parametric (1: spm) or non-parametric (2: SnPM) analysis. SnPM is
% limited to two-way interactions, so if three factors are select for
% interaction with SnPM it will output all 2-way interactions (actually, t-tests on subtracted data)

S.para = 2;
S.folder =0;
S.spmstats_path = S.data_path;
% factors and statistical model
S.factors = {'Grp', 'Subject'}; % must include a subject factor at the end
S.factortype = {'g','s'}; % w = within, s = subject, g = subject group
S.useIDfile = 0; % use data from Participant ID data file (excel)?
S.imglist = inputNifti;
S.cond_list = 1;
S.cov_names = {};
S.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
S.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA
% the following are for SnPM, not SPM
S.nPerm = 5000; % permutations
S.vFWHM = [20 20 20]; % variance smoothing (should be same as data smoothing used)
S.bVolm = 1; % 1=high memory usage, but faster
S.ST_U = 0.05; % cluster forming threshold
% Write residuals? For normality tests
S.resid = 0;
S.grpind = grpind;
S.interactions = [0 0]; % one column per factor; one row per interaction
% run design_batch function for each contrast
for iCon = nContrasts:-1:1,
    if sum(contrasts(iCon,:))>0 % no group contrast
        S.maineffects = [0 0]; % one column per factor 
        S.anapref = 'OneSampleT'; % directory prefix
        S=design_batch(S);
    elseif sum(contrasts(iCon,:))==0 % group contrast
        S.maineffects = [1 0]; % one column per factor 
        S.anapref = 'TwoSampleT'; % directory prefix
        S=design_batch(S);
    end
    %% Retrieve results
    pth = S.spm_path;
    
    niifile = {'lP+.img','lP-.img'};
    for ni = 1:length(niifile)
        V = load_nii(fullfile(pth,niifile{ni}));
        Y=V.img;
        Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
        Y=10.^(-Y); % for p values
        ptmp(:,iCon,ni) = Y;
    end
    niifile = {'lP_FWE+.img','lP_FWE-.img'};
    for ni = 1:length(niifile)
        V = load_nii(fullfile(pth,niifile{ni}));
        Y=V.img;
        Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
        Y=10.^(-Y); % for p values
        fweptmp(:,iCon,ni) = Y;
    end
    niifile = {'lP_FDR+.img','lP_FDR-.img'};
    for ni = 1:length(niifile)
        V = load_nii(fullfile(pth,niifile{ni}));
        Y=V.img;
        Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
        Y=10.^(-Y); % for p values
        fdrptmp(:,iCon,ni) = Y;
    end
    niifile = {'snpmT+.img','snpmT-.img'};
    for ni = 1:length(niifile)
        V = load_nii(fullfile(pth,niifile{ni}));
        Y=V.img;
        Y(~isfinite(Y))=[]; %delete NaN values from vector Y.
        Ttmp(:,iCon,ni) = Y;
    end
end%for
end%perform_glm_with_snpm
% [EOF]