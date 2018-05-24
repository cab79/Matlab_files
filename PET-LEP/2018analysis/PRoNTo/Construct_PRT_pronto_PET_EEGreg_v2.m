%% setup design batch for factorial design and estimation
clear all
addpath(genpath('C:\Data\Matlab\PRoNTo_dev-2.0.1'));
addpath(genpath('C:\Data\Matlab\SPM12'));
addpath('C:\Data\Matlab\Matlab_files');
%files required: participant data file with columns headed 'Subject' (should be characters, e.g. S1),
%'Group' (can be described by numbers or characters, but numbers recommended), 'Include' (must be numbers 
%with 0 meaning to exclude subject from analysis) and if you have covariates, columns headed with each
%covariate name. There should be just one header row, and one row for each
%subject.

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% name and location of the current design-batch file
D.batch_path = 'C:\Data\Matlab\Matlab_files\PET-LEP\2018analysis\PRoNTo\Construct_PRT_pronto_PET_EEGreg.m';
% template flexible factorial matlabbatch
D.batch = 'C:\Data\Matlab\Matlab_files\PET-LEP\2018analysis\PRoNTo\PRT.mat';
% root directory in which subject-specific folders are located
D.data_path = 'F:\Dell\bloodA\Image_analysis_files\examples';
% directory in which image masks are saved
D.mask_path = '';
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
D.pdatfile = 'C:\Data\PET-LEP\Participant_data.xlsx';
% names of headers in the above xls file:
    D.subhead = 'PET_ID';
    D.grphead = {'Group'}; 
    D.inchead = 'Include_PET';
% directory in which analyses will be saved (new folder created)
D.spmstats_path = 'C:\Data\PET-LEP\PET\pronto';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix and suffix of subject folder names (within 'data_path') either side of subject ID
D.anapref = ''; %directory prefix for this specific analysis
D.subdirpref = ''; % generic prefix
D.subdirsuff = '\PET'; % generic suffix
D.folder =1; % Is the data in a subject-specific folder?
D.identifier='_krr_GFP_pain-nopain_perm1000_log_log_new'; % optional identifer to add to end of outputted SPM folder name

% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% list of image names within each subject folder
D.imglist = {};
% alternaitvely, name of column headers in Participants_data file
% containing image names
D.imglist_columnheaders = {
            'pain_scan'
            'nopain_scan'
            }; 
D.imgpref = 'wr*';
D.imgsuff = '_nf32_i12_nnls_Vd_mag-310_rsl.img';

%% analysis design and parameters
%-------------------------------------------------------------
D.pronto = 1; % multivariate
% specify a time window to analyse
%D.time_ana = [-3000 -500]; % applies a mask to the data
D.time_ana = []; % applies a temporal mask to the data (first level)
D.maskfile = 'C:\Data\Matlab\mricron\templates\ch2better.nii'; % alternatively, directly specify a mask file (only if D.time_ana is empty)
%D.time_ana = [0 1500]; % applies a temporal mask to the data (first level)
%D.timewin = [];% apply windowing over the range of D.time_ana? Provide window size
D.ROIfile = 'C:\Data\Matlab\spm12\toolbox\aal\aal2.nii';

% contrast: regression will take place on a single image per subject. So,
% if there are multiple images (multiple conditions), new image will be
% calculated first based on a contrast. Use 1s and 2s, not -1 and 1.
% E.g. [1 1; 1 2] - first column states how many individual contrasts are
% performed, second column is the two images that are contrasted.
D.regress_contrast =  [
    1 1
    1 2
    ];

% For SnPM or Pronto, if interactions enter the design then subtractions
% are needed by specifying 'contrast' here
D.fileoptype = 'contrast';
%D.fileoptype = 'meancond';
% overwrite previous images with same name
D.overwrite =0;

D.regress_contrastname =  'name';
D.grp_list = [1 2];

%D.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
%D.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA

% names of nuisance covariates - NOT IMPLEMENTED FOR PRONTO YET
%cov_names = {'Age','Gender'};
D.cov_names = {
   % 'event1_peak1'
    'event1_peak2'
    %'event1_peak3'
    %'event1_peak4'
    %'event1_peak5'
    %'event1_peak6'
    %'event1_peak7'
    %'event1_peak8'
    %'event1_peak9'
    %'event1_peak10'
    %'event1_peak11'
    %'event1_peak12'
    %'event1_peak13'
    %'event1_peak14'
    %'event1_peak15'
   % 'event1_peak16'
   % 'event2_peak1'
    %'event2_peak2'
    %'event2_peak3'
    %'event2_peak4'
    %'event2_peak5'
   % 'event2_peak6'
   % 'event2_peak7'
    %'event2_peak8'
    %'event2_peak9'
    %'event2_peak10'
    %'event2_peak11'
    %'event2_peak12'
    %'event2_peak13'
    %'event2_peak14'
    %'event2_peak15'
    %'event2_peak16'
    };
% if there are many covariates, these will be included in one multiple
% regression, unless this next value is set to 1, in which case separate
% regressions will be run.
D.runmultiple = 1;

% the following are for spm analysis, not snpm
%D.GMsca = [0 0 0]; %grand mean scaling
%D.ancova = [0 0 0]; %covariate
% after model estimation, constrasts to display
%D.fcontrasts = {
%    };

% use kernel?
D.kernel = 1;

% permutation testing
D.permtest = 0;
D.saveallweights = 0; % Requires added code to prt_compute_weights_class.m, line 287:
                %if length(d.coeffs)>size(d.datamat,1)
                %    keep_idx = train_idx(find(ID(:,6)==1));
                %    d.coeffs = d.coeffs(train_idx);
                %end

% Data operations
% 1. Sample averaging (within blocks): constructs samples by computing the average of all
% scans within each block or event for each subject and scondition.
% 2. Sample averaging (within subjects): constructs samples by computing the average of
% all scans within all blocks for each subject and scondition.
% 3. Mean centre features using training data: subtract the voxel-wise mean from each
% data vector.
% 4. Divide data vectors by their norm: scales each data vector (i.e. each example) to lie
% on the unit hypersphere by dividing it by its Euclidean norm.
D.data_op  = {2 3}; 
% mean centring (same as 3 above)
%D.meancentre = 0;

D.transform = 'log';

% regression machine options:
%D.machine = 'gpr';
D.machine = 'krr';

% classification machine options:
%D.machine = 'svm_binary';
%D.machine = 'gpc_binary';
%D.machine = 'gpc_multi'; % PRT_MODEL.M HAS BEEN MODIFIED TO ONLY ALLOW BINARY
%D.machine = 'mkl';

%D.cv_type = 'cv_lkso';D.nfolds = 2;
D.cv_type = 'cv_loso';

%% run design_batch function
prtstats = struct;
if D.runmultiple
    D.cov_names_orig = D.cov_names;
    for c = 1:length(D.cov_names_orig)
        clear functions PRT
        D.cov_names = D.cov_names_orig(c);
        Dout=design_batch(D);
        load(fullfile(Dout.spm_path,'PRT.mat'));
        if c==1
            prtstats = PRT.model(1).output.stats;
        else
            prtstats(c) = PRT.model(1).output.stats;
        end
    end
else
    design_batch(D);
end

if isfield(prtstats,'permutation')
    for c = 1:length(D.cov_names_orig)
        prtstats(c).permutation_pval_r2 = prtstats(c).permutation.pval_r2;
    end
end

dt = datestr(datetime,30);
save(fullfile(D.spmstats_path,['prtstats_' dt]),'prtstats');

prt
