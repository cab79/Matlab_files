% List of open inputs
% Source inversion: M/EEG datasets - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/scratch/cb802/Matlab_files/CRPS_digits/SPM batch/inverse_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Source inversion: M/EEG datasets - cfg_files
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
