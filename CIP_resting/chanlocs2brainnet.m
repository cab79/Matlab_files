function [nodes] = chanlocs2brainnet(datapath,basename,nodecol,nodesize)

%if isempty(datapath)
%    datapath = 'W:\Data\CRPS_resting\EEG\';
%end

%if isempty(basename)
%    files = dir([datapath '*.set']);
%    basename = files(1).name;
%end

[num,txt,raw] = xlsread('ROIs.xlsx');
ROInames = txt(3:162,7);

X = num(:,1);
Y = num(:,2);
Z = num(:,3);

nodes = [];
nodes = [num2cell(X) num2cell(Y) num2cell(Z) num2cell(nodecol') num2cell(nodesize') ROInames];

