filepath = 'C:\Data\CORE\PilotData'; % MUST include the '\' on the end
savepath = 'C:\Data\CORE\Preprocessed';
addpath(genpath('M:\Matlab\Matlab_files\CORE'));
parts = [2];
combine=0; % don't combine for CORE - can be combined after epoching

subj = dir(fullfile(filepath,'test*'));
for s = 1:length(subj)
    dataimport(filepath,[subj(s).name '*'],parts,combine,savepath)
end