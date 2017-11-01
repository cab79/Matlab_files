%%PLOT TOPOGRAPHIC MAPS

%File location
filepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Preprocessed';
files = dir('*_cleaned.set');
files_ana = 1:length(files);

cd(filepath)

%Plot channel ERP's in scalp array
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); % extract parts of filename for use later, e.g. saving
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath); 
    EEG = eeg_checkset( EEG );
    figure;pop_topoplot(EEG,1, 299,'',[1 1] ,0,'electrodes','on', 'maplimits', [-0.7, 0.7]);
end
