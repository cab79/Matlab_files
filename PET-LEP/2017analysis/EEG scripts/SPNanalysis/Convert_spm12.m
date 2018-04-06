clear all
filepath = 'C:\Data\PET-LEP\Preprocessed';
cd(filepath)
files = dir('*cleaned.set');
eventtypes = {'S  1','S  2'}; use_etype = [1 2]; 
flip=0;

for f = 1:length(files)
    EEG = pop_loadset(files(f).name,filepath);
    [pth nme ext] = fileparts(EEG.filename);
    [pth orig_datfile ext] = fileparts(EEG.filename)
    orig_datfile = [orig_datfile '.fdt'];
    EEG.datfile = [nme '.fdt'];
    if isfield(EEG,'dataset'); EEG.dataset = [nme ext];end
    if exist('fullfile(filepath, orig_datfile)','file')
        movefile(fullfile(filepath, orig_datfile),fullfile(filepath, EEG.datfile));
    end
    [pth nme ext] = fileparts(EEG.datfile);
    EEG.dataset = [nme '.set'];
    prefix = 'spm12_';
    if flip==1 && (~isempty(strfind(nme, 'right')) || ~isempty(strfind(nme, 'Right'))); 
        EEG = flipchan(EEG); 
        prefix = 'spm12_flip_';
    end
    EEG.outfile = [prefix spm_str_manip(EEG.dataset,'tr')];
    
    fevents = zeros(1,size(EEG.data,3));
    for ev = use_etype
        fevents(1,find(strcmp({EEG.epoch.eventtype},eventtypes{use_etype(ev)}))) = use_etype(ev);
    end
    
    EEG.conditionlabels = cellfun(@num2str, num2cell(fevents), 'UniformOutput', false);
    
    EEG.mode ='epoched';
    EEG.timewin = [-3500 1998];
    spm_eeg_convert(EEG);
    %end
    clear EEG;
end
matlabmail
%x=1;save('x.mat','x');