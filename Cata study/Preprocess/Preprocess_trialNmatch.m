clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned.set');
load('C:\Data\Catastrophising study\Orig\chanlocs.mat');

stimtypematch = {'c1','c3';
                 'c2','c4';
                 'c5','c7';
                 'c6','c8'};

files_ana = [20]%1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    for c = 1:length(stimtypematch)
        
        EEGc = pop_selectevent(EEG, 'type', stimtypematch{c,1});
        EEGx = pop_selectevent(EEG, 'type', stimtypematch{c,2});
        EEGnc = pop_selectevent(EEG, 'type', stimtypematch{c,1},'invertepochs','on');
        if EEGc.trials>EEGx.trials
            EEGc = pop_select(EEGc, 'trial', randsample(1:size(EEGc.data,3), EEGx.trials));
        end
        EEG = pop_mergeset(EEGnc,EEGc);
    end
  
    sname = [C{1} '_' C{2} '_cleaned_trialNmatch.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
end