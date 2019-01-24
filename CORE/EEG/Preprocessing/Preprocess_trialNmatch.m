clear all
filepath = 'C:\Data\CORE\EEG\ana\prep\cleaned\part2';
cd(filepath);
files = dir('*2_merged_cleaned.set');
load('C:\Data\Matlab\Matlab_files\CORE\EEG\Supporting_functions\chanlocs.mat');


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


files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    ALLEEG=EEG;
    for c = 1:length(stimtypematch)
        % select event type
        selectepochs_c = find(conds==stimtypematch(c,2));
        selectepochs_x = find(conds==stimtypematch(c,1));
        
        %EEGnc = pop_select(EEG_ALL, 'trial', stimtypematch{c,1},'invertepochs','on');
        if length(selectepochs_c)>1 && length(selectepochs_x)>1
            EEGc = pop_select(EEG,'trial',selectepochs_c);
            EEGx = pop_select(EEG,'trial',selectepochs_x);
            if EEGc.trials>EEGx.trials
                EEGc = pop_select(EEGc, 'trial', randsample(1:size(EEGc.data,3), EEGx.trials));
            elseif EEGx.trials>EEGc.trials
                EEGx = pop_select(EEGx, 'trial', randsample(1:size(EEGx.data,3), EEGc.trials));
            end
            ALLEEG(c)=pop_mergeset(EEGx,EEGc);
        elseif ~isempty(selectepochs_c) || ~isempty(selectepochs_x)
            ALLEEG(c)=pop_select(EEG,'trial',sort([selectepochs_c selectepochs_x]));
        end
    end
    ALLEEG(cellfun(@isempty,{ALLEEG.setname}))=[];
    EEG = pop_mergeset(ALLEEG, 1:length(ALLEEG));
    sname = [C{1} '_' C{2} '_cleaned_tm.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
end