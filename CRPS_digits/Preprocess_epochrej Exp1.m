clear all
if isunix
    filepath = '/scratch/cb802/Data';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end
raw_path = fullfile(filepath,'CRPS_raw/Raw');
cd(raw_path);
files = dir('*Exp1*t.set');
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.
load chanlocs

for f = 1:length(files)
   
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',raw_path);
    EEG = pop_reref( EEG, []);
    EEG = pop_subcomp( EEG, [], 0);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    EEG = pop_autorej(EEG, 'nogui','on','threshold',100,'startprob',5);
    EEG = pop_saveset(EEG,'filename',EEG.filename,'filepath',raw_path); 
end