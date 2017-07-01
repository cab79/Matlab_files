clear all
if isunix
    filepath = '/scratch/cb802/Data/FMOA_resting';
else
    filepath = 'W:\Data\FMOA_resting';
end
cd(filepath);
files = dir(fullfile(filepath,'EEG','*_total_data_ICA_ca*'));
filesb = dir(fullfile(filepath,'EEG','*_block_info*'));
conntype = 'eeglab_dwpli';
mattype = 'max';
num = 1;
EEG = pop_loadset('filename','ExampleEEGLAB.set','filepath',filepath);

for f = 1:length(files)
    if ~isunix
        g = length(files)-f+1-(num-1);
    else
        g = f+num-1;
    end
    filename = files(g).name;
    filenameb = filesb(g).name;
    [pth basename ext] = fileparts(filename);
    [pth basenameb ext] = fileparts(filenameb);
    savename = fullfile(filepath,conntype,[basename '_' conntype '_randgraph.mat']);
    %savename = fullfile(filepath,conntype,[basename '_' conntype '_graph.mat']);

    if ~exist(savename,'file');
        %savepath = fullfile(filepath,conntype);
        %eeglabcoherence_FMOA(filepath,basename,basenameb,EEG,conntype)
        calcgraph_FMOA(basename,conntype,mattype,'randomise','on','latticise','off');
        %calcgraph_FMOA(basename,conntype,mattype,'randomise','off','latticise','off');
    end
end

matlabmail