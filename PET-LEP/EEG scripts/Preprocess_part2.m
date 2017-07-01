clear all
filepath = 'C:\Data\PET-LEP\Preprocessed';
cd(filepath);
files = dir('*ICA.set');
load('C:\Data\PET-LEP\Preprocessed\chanlocs.mat');

ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebin = [-3.5 -3];
stimtypes = {'S  1','S  2'};

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    EEG = pop_subcomp( EEG, [], 0);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    EEG = FTrejman(EEG,[0 0]);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    EEG = pop_reref( EEG, []);
    %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',15);
    sname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

    EEGall=EEG;
    if ALLEEG_save==1
        ALLEEG=struct;
    end
    for st = 1:length(stimtypes)
        EEG = pop_selectevent(EEGall,'type',stimtypes{st});
        if ALLEEG_save==1
            if st==1
                ALLEEG = EEG; 
            else 
                ALLEEG(st)=EEG;
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