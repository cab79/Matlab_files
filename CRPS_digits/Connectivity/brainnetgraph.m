%function brainnetgraph(datapath,basename,graphdir,freq,thresh,plot)
plot=1;
%if isempty(datapath)
    datapath = 'W:\Data\CRPS_resting\EEG\';
%end

%if isempty(basename)
    files = dir([datapath '*subcomp.set']);
    [pth basename ext] = fileparts(files(1).name);
%end

%if isempty(graphdir)
    graphdir = 'ft_dwpli';
%end

%if isempty(thresh)
    thresh = 0.5;
%end

%if isempty(freq)
    freq = 3;
%end

tvals = 0.5:-0.025:0.1;

ti = find(tvals==thresh);

load(fullfile(datapath,graphdir,[basename '_ftdwpli_graph.mat'])); 
load(fullfile(datapath,graphdir,[basename '_ftdwpli.mat'])); 

ftmat = squeeze(matrix(freq,:,:));
weimat = threshold_proportional(ftmat,thresh);
binmat = double(threshold_proportional(ftmat,thresh) ~= 0);

cen = graphdata{find(strcmp(graphdata(:,1),'centrality')),2};
cen = squeeze(cen(freq,ti,:));
nodesize = cen;

nodecol = graphdata{find(strcmp(graphdata(:,1),'modules')),2};
nodecol = squeeze(nodecol(freq,ti,:));



nodes = chanlocs2brainnet(datapath,basename,nodecol,nodesize);

fid = fopen(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '.node']),'w');
[nrows,ncols] = size(nodes);
for row = 1:nrows
    fprintf(fid,'%f   %f   %f   %d   %f   %s\n',nodes{row,1},nodes{row,2},nodes{row,3},nodes{row,4},nodes{row,5},nodes{row,6});
end
fclose(fid);

fid = fopen(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_weighted.edge']),'w');
[nrows,ncols] = size(nodes);
fprintf(fid,[repmat('%g ',1,size(weimat,2)),'\n'], weimat);
fclose(fid);

fid = fopen(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_binary.edge']),'w');
[nrows,ncols] = size(nodes);
fprintf(fid,[repmat('%g ',1,size(binmat,2)),'\n'], binmat);
fclose(fid);

%save(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_weighted.edge']),'weimat','-ascii');
%save(fullfile(datapath,graphdir,[basename '_' graphdir '_graph_freq' num2str(freq) '_thresh_' num2str(ti) '_binary.edge']),'binmat','-ascii');

if (plot == 1) BrainNet_MapCfg('W:\Data\CRPS_resting\EEG\ft_dwpli\H07.100Hz.Exp3_BELB_subcomp_ft_dwpli_graph_freq3_thresh_1_binary.edge',...
        'W:\Data\CRPS_resting\EEG\ft_dwpli\H07.100Hz.Exp3_BELB_subcomp_ft_dwpli_graph_freq3_thresh_1.node',...
        'W:\Data\CRPS_resting\EEG\Cfg.mat'); end
