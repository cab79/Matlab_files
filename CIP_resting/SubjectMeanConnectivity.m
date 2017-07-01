clear all
if isunix
    filepath = '/scratch/cb802/Data/CIP_resting/fMRI';
else
    filepath = 'W:\Data\CIP_resting\fMRI';
end
cd(filepath);
conntype = 'fMRI_pearsons';
mattype = [];
anatype = '';
cd(fullfile(filepath,conntype,anatype));
files = dir(fullfile(filepath,conntype,anatype,['*' conntype '.mat']));
healthy = 1:39;
patients = 40;
meanfiles = 1:39;
filenamesave = 'Healthy_';
filetypes = {''};

matrix = [];
bootmat = [];
for f = meanfiles
    for ft = 1:length(filetypes)
        filename = files(f).name;
        [pth basename ext] = fileparts(filename);
        datfile = dir(fullfile(filepath,conntype,anatype,[basename '.mat']));
        dat = load(datfile.name);
        matrix = cat(3,dat.matrix, matrix);
        bootmat = cat(3,dat.bootmat, bootmat);
    end
end

matrix = squeeze(mean(matrix,3));
bootmat = squeeze(mean(bootmat,3));

save(fullfile(filepath, conntype, [filenamesave conntype '.mat']),'matrix','bootmat');

calcgraphRho(filenamesave,conntype,mattype,'randomise','on','latticise','off');
calcgraphRho(filenamesave,conntype,mattype,'randomise','off','latticise','off');