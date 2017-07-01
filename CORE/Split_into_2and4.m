clear all
filepath = 'C:\Data\CORE\Preprocessed';
cd(filepath);
files = dir('*orig.set');
ISIs = [1, 0.4];
addpath(genpath('M:\Matlab\Matlab_files\CORE\Supporting functions'));
files_ana = [6]% 1:length(files);
markers = {'STIM' 'DIN2' 'break cnt'};
use_marker = 3;
mult=1.1;
partnames = {'2' '4'};

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    if strcmp(C{2},'24')
        EEG = pop_loadset('filename',orig_file,'filepath',filepath);
        markidx = strcmp({EEG.event.type},markers(use_marker)); 
        
        if strcmp(markers(use_marker),'break cnt')
            marki = find(markidx);
            if length(marki)==2
                latbreak = EEG.event(marki(2)).init_time;
            else
                error('too many breakcnt events');
            end
        else
            lats = [EEG.event(markidx).latency];
            latdiff = lats(2:end)-lats(1:end-1);
            breaks = [];
            for i = 2:length(latdiff)-1
                if latdiff(i-1)<1000*ISIs(1)*mult && latdiff(i-1)>1000*ISIs(1)*(2-mult)
                    if latdiff(i)>1000*ISIs(1)*mult
                        if latdiff(i+1)<1000*ISIs(2)*mult && latdiff(i+1)>1000*ISIs(2)*(2-mult)
                            breaks = [breaks i];
                        end
                    end
                end
            end

            if length(breaks)==1
                mark_event = EEG.event(markidx);
                latbreak = mark_event(breaks).init_time+2;
            else
                error('too many breaks');
            end
        end
        
        EEGall = EEG;
        
        all_lat = [0 latbreak EEGall.event(end).init_time+2];
        for parts = 1:2
            timrange = [all_lat(parts) all_lat(parts+1)];
            EEG = pop_select(EEGall, 'time', timrange);
            EEG.filename = [C{1} '_' partnames{parts} '_orig.set'];
            EEG.setname = EEG.filename;
            EEG = pop_saveset(EEG,'filename',EEG.filename,'filepath',filepath); 
            clear EEG
        end
        
    end
end