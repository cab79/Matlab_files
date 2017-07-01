clear all
filepath = 'C:\Data\CHEPS\Preprocessed';
cd(filepath);
files = dir('*_1st_ICA.set');
load('M:\Matlab\Matlab_files\CORE\Supporting functions\chanlocs.mat');

ALLEEG=struct;
ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebins = [-0.5 0]; % for epoching, TSOT(4)
        
files_ana = [1]% 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    EEG = pop_reref( EEG, []);
    %EEG = pop_subcomp( EEG, [], 0);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    EEG = FTrejman(EEG,[0 0]);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',15);
    sname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
    
end

for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',[C{1} '_' C{2} '_cleaned.set'],'filepath',filepath);
    
    %if strcmp(C{2},'2')
        basebin = basebins(1,:);
    %elseif strcmp(C{2},'4')
    %    basebin = basebins(2,:);
    %end
    
    EEG = pop_rmbase(EEG, basebin*1000);
    
    ALLEEG=EEG;
    sname = [C{1} '_' C{2} 'ALLEEG.mat'];
    save(sname,'ALLEEG');
   
end