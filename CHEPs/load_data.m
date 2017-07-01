filepath = 'C:\Data\CHEPs\Raw'; % MUST include the '\' on the end
savepath = 'C:\Data\CHEPs\Preprocessed';
parts = [1];
combine=0; % don't combine for CORE - can be combined after epoching

subj = dir(fullfile(filepath,'CHEPS*'));
for s = 1:length(subj)
    dataimport(fullfile(filepath),[subj(s).name '*'],parts,combine,savepath)
end