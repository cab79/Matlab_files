%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_condor_fit_step1

pth=pwd;

% cleanup from last run
delete('input*.mat');
delete(fullfile(pwd,'Data/input_files','input*.mat'));
delete('output*.mat');
delete('*.out');
delete('*.log');
delete('*.err');

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'part2'; %
S.path.seq = fullfile(pth, 'Data'); % unprocessed data in original format
S.path.raw = fullfile(pth, 'Data'); % unprocessed data in original format
S.path.prep = fullfile(pth, 'Data'); % folder to save processed data
S.path.hgf = fullfile(pth, 'Data'); % folder to save processed data
S.path.design = fullfile(pth, 'Data'); % 
S.fname.parts = {'prefix','subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.groups = {};
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = fullfile(pth, 'Data','Participant_data.xlsx'); % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% unique save name extension
sname = datestr(now,30);
S.sname=sname;
if ~isempty(S.select.groups)
    S.sname = [S.sname '_' S.select.groups{:}];
end

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S); 

% response data processing
% NOW OCCURS AFTER MODEL SPECIFICATION TO ALLOW MULTIPLE ALPHAS
S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 1;
S.RT.min = 0.2; % min RT to consider
S.RT.max = 1;
S.RT.fill_non_response = 0; % RT value to add to trials in which there was no response: acts to include inaccurary info within RT data
S.save.tables = 0;
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% split data into training and testing sets (if we want to test for prediction of behaviour)
S.frac_train = 1; % set to 0 to include all data in training set AND test set
D_fit=D_prep;
if S.frac_train>0
    for d=1:length(D_prep)
        cond = D_prep(d).dt.design(2,:); % get conditions
        ucond = unique(cond);
        % random indices from each cond
        S.train_idx = [];
        for u = 1:length(ucond)
            cond_idx = find(cond==ucond(u));
            S.train_idx = [S.train_idx randsample(cond_idx,S.frac_train*length(cond_idx))];
        end
        S.test_idx = 1:length(D_fit(d).HGF.u);
        S.test_idx(S.train_idx)=[];
        %D_fit(d).HGF.u = D_prep(d).HGF.u(sort(S.train_idx),:);
        D_fit(d).HGF.y(S.test_idx) = nan;
    end
end

S.HGF.plottraj = 0; % turn off if doing multiple simulations!
S.parallel = 0; % parallelise over subjects

[nsub,nmods] = make_input_files(D_fit,S,pth); 

create_job_submission_file(pth,nsub,nmods)
    
% save info
fid = fopen(fullfile(pth,'nsub.txt'),'w');
fprintf(fid,'%d',nsub);
fclose(fid);
fid = fopen(fullfile(pth,'nmods.txt'),'w');
fprintf(fid,'%d\n',nmods);
fclose(fid);

quit 

function [nsub,nmods] = make_input_files(D,S,pth)

cd(pth)

if ~exist(fullfile(S.path.hgf,'input_files'),'dir')
    mkdir(fullfile(S.path.hgf,'input_files'))
end

% 1 - Model initialisation using pre-defined priors (ideally, Bayes optimal priors with wide variances)
S.perc_models=[2]; % 1 3 9 10 11 12
S.resp_models = [1:6]; % 2 12 21
nmods = [length(S.perc_models),length(S.resp_models)];
nsub = length(D);
ninputs = prod(nmods) * nsub;

for pm=1:length(S.perc_models)
    for rm=1:length(S.resp_models)
        S.perc_model = S.perc_models(pm); 
        S.resp_model = S.resp_models(rm); 
        S.prc_config = ['GBM_config_CORE_percmodel' num2str(S.perc_model) '_BOpriors']; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
        S=CORE_perceptual_models(S);
        S=CORE_response_models(S);

        prc_model = S.prc_config;
        pmodel = strrep(prc_model,'_config',''); 
        obs_model = S.obs_config;
        opt_algo = 'tapas_quasinewton_optim_config';

        S.use_y_col = find([any(strcmp(S.resp_modelspec.responses,'Ch')), any(strcmp(S.resp_modelspec.responses,'RT')),any(strcmp(S.resp_modelspec.responses,'EEG'))]);

        for d = 1:length(D)

            if ~isfield(S,'nstim')
                nst=length(D(d).HGF(1).u);
            else
                if isempty(S.nstim)
                    nst=length(D(d).HGF(1).u);
                else
                    nst=S.nstim;
                end
            end
            y = D(d).HGF(1).y(1:nst,:); 
            u = D(d).HGF(1).u(1:nst,:);

            ii = (pm-1)*length(S.resp_models)*length(D) + (rm-1)*length(D) + d;
            disp(['creating input file ' num2str(ii) '/' num2str(ninputs)])
            save(fullfile(S.path.hgf,'input_files',['input' num2str(ii-1) '.mat']),'y','u','prc_model','obs_model','opt_algo','S')

        end
    end
end

disp('FINISHED');
    
function create_job_submission_file(pth,ns,nmods)

disp('creating job submission file')
A = {
    'executable=CORE_condor_fit_job.exe'
    'indexed_input_files=input.mat'
    'indexed_output_files=output.mat'
    'indexed_stdout=CORE_condor_fit_job.out'
    'indexed_stderr=CORE_condor_fit_job.err'
    'indexed_log=CORE_condor_fit_job.log'
    'max_run_time=240'
    ['total_jobs=' num2str(ns*prod(nmods))]
};

fid = fopen(fullfile(pth, ['CORE_condor_fit_job_run.sub']),'w');
for i = 1:size(A,1)
    fprintf(fid,'%s\n',A{i});
end
fclose(fid);