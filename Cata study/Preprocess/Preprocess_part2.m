clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*_orig_ICA.set');
load('C:\Data\Catastrophising study\Orig\chanlocs.mat');

ALLEEG_save = 0; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebin = [-5.5 -5];
endbin = [1.5 2];
stimtypes = {'c0','c1','c2','c3','c4','c5','c6','c7','c8'};

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    EEG = pop_subcomp( EEG, [], 0); % remove selected ICA components
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    %EEG = pop_reref( EEG, []);
    EEG = FTrejman(EEG,[0 0]);
    
    % NEW DETREND
    %for i = 1:EEG.trials
    %    EEG.data(:,:,i) = detrend4(EEG.data(:,:,i)',1,250,250)';
    %end
   
    sname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

    if ALLEEG_save
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
end
if ALLEEG_save==2
    sname = ['Allfiles_nochange_ALLEEG.mat'];
    save(sname,'ALLEEG');
end