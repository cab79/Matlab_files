clear all
if isunix
    filepath = '/scratch/cb802/Data/CIP_resting/fMRI/RS_dosenbach_TimeSeries';
    savepath = '/scratch/cb802/Data/CIP_resting/fMRI';
else
    filepath = 'W:\Data\CIP_resting\fMRI\RS_dosenbach_TimeSeries';
    savepath = 'W:\Data\CIP_resting\fMRI';
end
cd(filepath);
conntype = 'fMRI_pearsons';
grplist = [1 2];
run('W:\Matlab_files\CIP_resting\loadsubj.m')
subjects = subjlists(grplist);
mattype = [];
trunc = 157;

for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        basename = subj;
        %savename = fullfile(savepath,conntype,[basename conntype '_randgraph.mat']);
        savename = fullfile(savepath,conntype,[basename conntype '_graph.mat']);

        if ~exist(savename,'file');
            %Create_pearsonsR_matrices(filepath,basename,conntype, trunc);
            calcgraphRho(basename,conntype,mattype,'randomise','off','latticise','off');
            %calcgraphRho(basename,conntype,mattype,'randomise','on','latticise','off');
        end
    end
end

matlabmail