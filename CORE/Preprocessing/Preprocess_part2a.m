%% NEW VERSION

clear all
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);

fname_ext = '_1st_ICA';
fname_ext2 = '';
%fname_ext2 = '_ACSTP';
files = dir(['*t' fname_ext fname_ext2 '.set']);
load('C:\Data\Matlab\Matlab_files\CORE\Supporting_functions\chanlocs.mat');

ALLEEG=struct;
ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebins = [-0.2 0; % for epoching, TSOT(2)
            -0.05 0]; % for epoching, TSOT(4)
        
files_ana = [87 88]% 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    EEG = pop_reref( EEG, []);
    EEG = pop_subcomp( EEG, [], 0);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    sname = [C{1} '_' C{2} fname_ext '_' C{3} '_cleaned' fname_ext2 '.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
end
    
complete = '';
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    if ~strcmp(complete,C{1}) 
        leftname = [C{1} '_' C{2} fname_ext '_left_cleaned' fname_ext2 '.set'];
        rightname = [C{1} '_' C{2} fname_ext '_right_cleaned' fname_ext2 '.set'];
        LEEG = pop_loadset('filename',leftname,'filepath',filepath);
        REEG = pop_loadset('filename',rightname,'filepath',filepath);
        EEG = pop_mergeset(LEEG,REEG);
        sname = [C{1} '_' C{2} fname_ext '_merged' fname_ext2 '.set'];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
        complete = C{1};
    end
end