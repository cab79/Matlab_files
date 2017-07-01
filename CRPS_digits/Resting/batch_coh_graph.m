clear all
if isunix
    filepath = '/scratch/cb802/Data/CRPS_resting/EEG';
else
    filepath = 'W:\Data\CRPS_resting\EEG';
end
cd(filepath);
files = dir(fullfile(filepath,'*100Hz.Exp3*subcomp.set'));
conntype = 'ftdwpli';
%conntypesave = 'ftdwpli';
mattype = '';
%nu = 26;
num = 15;

for f = 1:length(files)
    if isunix
        g = length(files)-f+1-(num-1);
    else
        g = f+num-1;
    end
    filename = files(g).name;
    [pth basename ext] = fileparts(filename);
    savename = fullfile(filepath,conntype,[basename '_' conntype '_randgraph.mat']);
    %savename = fullfile(filepath,conntype,[basename '_' conntype '_graph.mat']);
    %savename = fullfile(filepath,conntype,[basename '_' conntype '.mat']);

    if ~exist(savename,'file');
        %ftcoherence(basename)
        calcgraph(basename,conntype,mattype,'randomise','edgerand','latticise','off');
        %calcgraph(basename,conntype,mattype,'randomise','on','latticise','off');
        %calcgraph(basename,conntype,mattype,'randomise','off','latticise','off');
    end
end

matlabmail