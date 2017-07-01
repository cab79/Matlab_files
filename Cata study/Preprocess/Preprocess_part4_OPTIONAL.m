clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned_ICA.set');
load('C:\Data\Catastrophising study\Preprocessed\chanlocs.mat');

ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebin = [-3.5 -3];
stimtypes = {'c0','c1','c2','c3','c4','c5','c6','c7','c8'};

%files_ana = [39];%usable if further cleaned
%files_ana = [12,17,23,28,35];%warnings or unusable
%files_ana = [2:11 13:16 18 20 24 26 29 31:34 36:37 40]; % the rest
files_ana = 25%[1 19 38]; % no FTrejman

for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    EEG = pop_subcomp( EEG, [], 0); % remove selected ICA components
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    %EEG = FTrejman(EEG,[0 0]);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    EEG = pop_reref( EEG, []);
    %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',15);
    sname = [C{1} '_' C{2} '_cleaned2.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

    EEGall=EEG;
    if ALLEEG_save==1
        ALLEEG=struct;
    end
    for st = 1:length(stimtypes)
        if any(strcmp({EEGall.event.type},stimtypes{st}))
            EEG = pop_selectevent(EEGall,'type',stimtypes{st});
            if ALLEEG_save==1
                if st==1
                    ALLEEG = EEG; 
                else 
                    ALLEEG(st)=EEG;
                end
            end
        end
    end
    if ALLEEG_save==1
        sname = [C{1} '_' C{2} '_ALLEEG.mat'];
        save(sname,'ALLEEG');
    end
    
    if ALLEEG_save==2
        if f==files_ana(1)
            ALLEEG = EEG; 
        else 
            ALLEEG(length(ALLEEG)+1)=EEG;
        end
    end
end
if ALLEEG_save==2
    sname = ['Allfiles_nochange_ALLEEG.mat'];
    save(sname,'ALLEEG');
end