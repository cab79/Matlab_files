clear all
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);
files = dir('*4_merged_cleaned.set');
load('C:\Data\Matlab\Matlab_files\CORE\Supporting_functions\chanlocs.mat');


stimtypematch = [
    1,3;
    2,4;
    5,7;
    6,8;
    9,11; 
    10,12; 
    13,15;
    14,16;
    17,19;
    18,20;
    21,23;
    22,24];


files_ana = 35:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    ALLEEG=EEG;
    for c = 1:length(stimtypematch)
        % select event type
        selectepochs = find(conds==stimtypematch(c,2));
        EEGc = pop_select(EEG,'trial',selectepochs);
        selectepochs = find(conds==stimtypematch(c,1));
        EEGx = pop_select(EEG,'trial',selectepochs);
        
        %EEGnc = pop_select(EEG_ALL, 'trial', stimtypematch{c,1},'invertepochs','on');
        if EEGc.trials>EEGx.trials
            EEGc = pop_select(EEGc, 'trial', randsample(1:size(EEGc.data,3), EEGx.trials));
        end
        ALLEEG(c)=pop_mergeset(EEGx,EEGc);
    end
    EEG = pop_mergeset(ALLEEG, 1:length(ALLEEG));
    sname = [C{1} '_' C{2} '_cleaned_tm.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
end