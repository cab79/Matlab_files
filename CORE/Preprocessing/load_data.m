%clear all
filepath = 'C:\Data\CORE\raw'; % MUST include the '\' on the end
savepath = 'C:\Data\CORE\Orig\';
parts = [2];
combine=0; % don't combine for CORE - can be combined after epoching
%firstsample = 1;
%lastsample = 5214475;

subj = dir(fullfile(filepath,'CORE*'));
for s =[24]%:length(subj)
    dataimport(fullfile(filepath,subj(s).name),[subj(s).name '*'],parts,combine,savepath);%,firstsample,lastsample)
end