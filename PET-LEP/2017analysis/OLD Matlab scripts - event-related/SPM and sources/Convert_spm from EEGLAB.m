clear all
if isunix
    loadpaths_unix
else
    loadpaths
end
cd(filepath)
files = dir('H07.Exp1.Right.set');

for f = 1:length(files)
    EEG = pop_loadset(files(f).name,filepath);
    
    EEG.filename = strrep(EEG.filename, '.left', '_left');
    EEG.filename = strrep(EEG.filename, '.Left', '_left');
    EEG.filename = strrep(EEG.filename, '.right', '_right');
    EEG.filename = strrep(EEG.filename, '.Right', '_right');
    EEG.filename = strrep(EEG.filename, '.flip', '_flip');
    EEG.filename = strrep(EEG.filename, '.Flip', '_flip');
    EEG.filename = strrep(EEG.filename, '.aff', '_aff');
    EEG.filename = strrep(EEG.filename, '.Aff', '_aff');
    EEG.filename = strrep(EEG.filename, '.Unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '.unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '_Left', '_left');
    EEG.filename = strrep(EEG.filename, '_Right', '_right');
    EEG.filename = strrep(EEG.filename, '_Flip', '_flip');
    EEG.filename = strrep(EEG.filename, '_Aff', '_aff');
    EEG.filename = strrep(EEG.filename, '_Unaff', '_unaff');
    EEG.filename = strrep(EEG.filename, '.Exp1', '_Exp1');
    
    [pth nme ext] = fileparts(EEG.filename); 
    
    [pth orig_datfile ext] = fileparts(EEG.filename)
    orig_datfile = [orig_datfile '.fdt'];
    EEG.datfile = [nme '.fdt'];
    if isfield(EEG,'dataset'); EEG.dataset = [nme ext];end
    
    if exist('fullfile(filepath, orig_datfile)','file');
        movefile(fullfile(filepath, orig_datfile),fullfile(filepath, EEG.datfile));
    end
    
    spm_eeg_convert_ui(EEG,'justread');
    clear EEG;
end