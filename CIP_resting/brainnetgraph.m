%function brainnetgraph(datapath,basename,graphdir,freq,thresh,plot)
plot=0;
%if isempty(datapath)
    datapath = 'W:\Data\CIP_resting\fMRI\fMRI_pearsons\';
%end

%if isempty(basename)
    files = dir([datapath '*pearsons.mat']);
    [pth basename ext] = fileparts(files(40).name);
%end

%if isempty(graphdir)
    graphdir = '';
%end

%if isempty(thresh)
    thresh = 0.1;
%end

%if isempty(freq)
    %freq = 3;
%end

%tvals = 0.5:-0.025:0.1;

% nodes of interest - save separate files for modules they are in
nodROI = [70 76]; % L and R post insula
nodnme = 'post_ins';
%nodROI = [26 28]; % L and R ant insula
%nodROI = [85 87 94 105 112 132]; % L and R precuneus

load(fullfile(datapath,graphdir,[basename '_graph.mat'])); 
load(fullfile(datapath,graphdir,[basename '.mat'])); 

ti = find(tvals==thresh);

ftmat = squeeze(matrix(:,:));
weimat = threshold_proportional(ftmat,thresh);
binmat = double(threshold_proportional(ftmat,thresh) ~= 0);

cen = graphdata{find(strcmp(graphdata(:,1),'centrality')),2};
cen = squeeze(cen(ti,:));
nodesize = cen;

nodecol = graphdata{find(strcmp(graphdata(:,1),'modules')),2};
nodecol = squeeze(nodecol(ti,:));

nodes = chanlocs2brainnet(datapath,basename,nodecol,nodesize);

nodenum = [];
for n = 1:length(nodROI)
    nodenum = [nodenum nodecol(nodROI(n))];
end
nodenum = unique(nodenum);

nodeidx = [];
for n = 1:length(nodenum)
    idx = find(nodecol==nodenum(n));
    nodeidx = [nodeidx idx];
end
nodeidx = sort(unique(nodeidx),'ascend');

nodes = nodes(nodeidx,:);
weimat = weimat(nodeidx,nodeidx);
binmat = binmat(nodeidx,nodeidx);

fid = fopen(fullfile(datapath,graphdir,[basename '_graph_thresh_' num2str(ti) '_' nodnme '.node']),'wt');
[nrows,ncols] = size(nodes);
for row = 1:nrows
    fprintf(fid,'%f   %f   %f   %d   %f   %s\n',nodes{row,1},nodes{row,2},nodes{row,3},nodes{row,4},nodes{row,5},nodes{row,6});
end
fclose(fid);

fid = fopen(fullfile(datapath,graphdir,[basename '_graph_thresh_' num2str(ti) '_' nodnme '_weighted.edge']),'wt');
[nrows,ncols] = size(nodes);
fprintf(fid,[repmat('%g\t',1,size(weimat,2)),'\n'], weimat);
fclose(fid);

fid = fopen(fullfile(datapath,graphdir,[basename '_graph_thresh_' num2str(ti) '_' nodnme '_binary.edge']),'wt');
[nrows,ncols] = size(nodes);
fprintf(fid,[repmat('%g\t',1,size(binmat,2)),'\n'], binmat);
fclose(fid);

%save(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_weighted.edge']),'weimat','-ascii');
%save(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_binary.edge']),'binmat','-ascii');

if (plot == 1) BrainNet_MapCfg('W:\Data\CIP_resting\fMRI\fMRI_pearsons\CON%00.S01.fMRI_pearsons_graph_thresh_1_binary.edge',...
        'W:\Data\CIP_resting\fMRI\fMRI_pearsons\CON%00.S01.fMRI_pearsons_graph_thresh_1.node',...
        'W:\Data\CIP_resting\fMRI\fMRI_pearsons\Cfg.mat'); end
