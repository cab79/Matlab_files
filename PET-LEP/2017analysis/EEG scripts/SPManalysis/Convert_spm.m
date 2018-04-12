clear all
filepath = 'C:\Data\PET-LEP\Preprocessed';
cd(filepath)
files = dir('*cleaned.set');
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
    if (~isempty(strfind(nme, 'left')) || ~isempty(strfind(nme, 'Left'))); 
        continue
    end
    %prefix = 'spm8_dc_';
    prefix = 'spm8_';
    if flip==1 && (~isempty(strfind(nme, 'right')) || ~isempty(strfind(nme, 'Right'))); 
        EEG = flipchan(EEG); 
        prefix = 'spm8_flip_';
        %prefix = 'spm8_flip_';
    end
    
    EEG.outfile = [prefix spm_str_manip(EEG.dataset,'tr')];
    spm_eeg_convert_ui(EEG,'justread');
    %end
    clear EEG;
end
%matlabmail
%x=1;save('x.mat','x');