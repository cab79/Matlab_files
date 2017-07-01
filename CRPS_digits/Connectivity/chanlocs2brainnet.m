function [nodes] = chanlocs2brainnet(datapath,basename,nodecol,nodesize)

if isempty(datapath)
    datapath = 'W:\Data\CRPS_resting\EEG\';
end

if isempty(basename)
    files = dir([datapath '*.set']);
    basename = files(1).name;
end

EEG = pop_loadset('filename',[basename '.set'],'filepath',datapath);
chanlocs = EEG.chanlocs;

nodes = [{chanlocs.X}' {chanlocs.Y}' {chanlocs.Z}' num2cell(nodecol) num2cell(nodesize) {chanlocs.labels}'];

