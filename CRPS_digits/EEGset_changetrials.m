clear all
close all

%---settings---%
%grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';%Exp2
%grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\';%Exp1 left v right
grplist = [51 52 53 54]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\';%Exp1 left v righ
eventtypes = {'FNUM','CNUM','ACCU'; [1:5], [0 1], [0 1]}; use_etype = [1 2];
allcond = [2 4 6 8 10];
selectcond = {2 [4 6 8] 10};
condnme = {'D1','M3','D5'};
run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath)
subjects = subjlists(grplist);
flip=1;

% create DATs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        if flip==1 && ~isempty(strfind(subj,'right'))
            EEG = flipchan(EEG); 
        end
        chanlocs=EEG.chanlocs;
        dsize = size(EEG.data);
        %EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        condevents = EEG2condevents(EEG,'STIM', eventtypes,use_etype);
        selectevents = ismember(condevents,allcond);
        EEGall = pop_select(EEG,'trial',find(selectevents));
        nevents=condevents(selectevents);
        sevents = unique(nevents);
        for se = 1:length(selectcond)
            condind = ismember(nevents,selectcond{se});
            EEG = pop_select(EEGall,'trial',find(condind));
            EEG = pop_saveset(EEG,'filename',[subj '_' condnme{se} '_change.set'],'filepath',fullfile(filepath,'LW'));
        end
    end
end


