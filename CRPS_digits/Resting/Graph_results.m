clear all
if isunix
    filepath = '/scratch/cb802/Data/CRPS_resting/EEG';
else
    filepath = 'W:\Data\CRPS_resting\EEG';
end
cd(filepath);
conntype = 'ftdwpli';
anatype = 'Edgerand';
cd(fullfile(filepath,conntype,anatype));
files = dir(fullfile(filepath,conntype,anatype,'*BELB*_graph.mat'));
healthy = 1:13;
patients = 14:25;
filetypes = {'BELB','LHAN','RELX','RHAN'};
nfreq = 5;
ncomp = 3;
ngraph = 9;
%ele = [80 79 78 77 76 61 66 70 73 72 60 65 69 59 64 68 58 63]; % wide definition (C, CP, P, PO)
%ele = [61 66 70 73 72 60 65 69 59 64 68 58 63]; % including PO (CP, P, PO)
%ele = [80 79 78 77 76 61 66 70 73 60 65 69 59 64]; % including central (C, CP, P)
%ele = [61 66 70 73 60 65 69 59 64]; % narrow definition (CP and P)
%ele = [60 65 69 59 64]; % very narrow definition (P only)
%ele = [72 60 65 69 59 64 68 58 63]; % (P,PO)
%ele = [80 79 78 77 76 61 66 70 73]; % including central (C, CP)
%ele = [61 66 70 73]; % CP only
ele = [61 66 70 73 60 69 59 64]; % narrow definition (CP and P) 66 0.008, 60 0.009, 65 0.003, 64 0.008
nonele = 1:92;
nonele(ele) = [];
aff_hand = [2 1 1 2 1 1 2 1 1 1 1 1 1 2 1 1 2 1 1 2 1 1 1 1 1];
aff_inx = [];
for af = 1:length(aff_hand)
    if aff_hand(af)==1
        affi = [1 2 3 4] + (af-1)*4;
    else
        affi = [1 4 3 2] + (af-1)*4;
    end
    aff_inx = cat(2, aff_inx, affi);
end

clear graphresults
graphresults{1,1} = 'clustering';
graphresults{1,2} = 'characteristic path length';
graphresults{1,3} = 'global efficiency';
graphresults{1,4} = 'modularity';
graphresults{1,5} = 'modules';
graphresults{1,6} = 'centrality';
graphresults{1,7} = 'modular span';
graphresults{1,8} = 'participation coefficient';
graphresults{1,9} = 'connection density';

graphresults = reshape(repmat(graphresults',1,nfreq*ncomp)',1,nfreq*ncomp*ngraph);
subhead1{1,1} = 'delta';
subhead1{1,2} = 'theta';
subhead1{1,3} = 'alpha';
subhead1{1,4} = 'beta';
subhead1{1,5} = 'gamma';
subhead1 = repmat(reshape(repmat(subhead1',1,ncomp)',1,nfreq*ncomp),1,ngraph);
graphresults(2,:) = subhead1;

subhead2{1,1} = 'allnodes';
subhead2{1,2} = 'r_par';
subhead2{1,3} = 'r_par_ratio';
subhead2 = repmat(subhead2,1,nfreq*ngraph);
graphresults(3,:) = subhead2;

graphresults(:,2:(nfreq*ncomp*ngraph+1)) = graphresults(:,1:(nfreq*ncomp*ngraph));
graphresults(1:3,1) = cell(3,1);

for f = 1:length(files)
    for ft = 1:length(filetypes)
        Rindstart = (f-1)*length(filetypes) + (ft-1) + 1;
        filename = files(f).name;
        [pth basename ext] = fileparts(filename);
        graphresults{3+Rindstart,1} = [basename(1:3) '_' filetypes{ft}];
        datfile = dir(fullfile(filepath,conntype,anatype,[basename(1:3) '*' filetypes{ft} '*_graph.mat']));
        randfile = dir(fullfile(filepath,conntype,anatype,[basename(1:3) '*' filetypes{ft} '*_randgraph.mat']));
        dat = load(datfile.name);
        rand = load(randfile.name);
        for gt = 1:9
            clear datamat
            Cindstart = nfreq*ncomp*(gt-1)+1;
            gdat = dat.graphdata{gt,2};
            rdat = rand.graphdata{gt,2};
            %if ~isreal(gdat) 
            %    gdat = abs(gdat);
            %    rdat = abs(rdat);
            %end
            rdat = mean(rdat,length(size(rdat)));
            rdat = squeeze(mean(rdat,2));
            gdat = squeeze(mean(gdat,2));
            normdat = gdat./rdat; 
            %normdat = normdat(isfinite(normdat));
            %fin = isfinite(normdat);
            %in_ind = find(fin==1);
            if size(normdat,2)>1
                
                datamat(1,:) = mean(normdat(:,nonele),2)';
                datamat(2,:) = mean(normdat(:,ele),2)';
                datamat(3,:) = datamat(2,:)./datamat(1,:);
            else 
                datamat = [normdat';nan(1,nfreq);nan(1,nfreq)];
            end
            datamat = reshape(datamat,1,nfreq*ncomp);
            resultsmat(Rindstart,Cindstart:Cindstart+(ncomp*nfreq-1)) = datamat;
        end
    end
end

resultsmat = resultsmat(aff_inx,:); % make LHAN attention to affected hand

graphresults(4:3+Rindstart,2:nfreq*ncomp*ngraph+1) = num2cell(resultsmat);

sizer = size(resultsmat);
resh = reshape(resultsmat,4,sizer(1)*sizer(2)/4);
resultsmat_mean = reshape(nanmean(resh,1),sizer(1)/4,sizer(2));

ft1 = reshape(resh(1,:),sizer(1)/4,sizer(2));
ft2 = reshape(resh(2,:),sizer(1)/4,sizer(2));
ft3 = reshape(resh(3,:),sizer(1)/4,sizer(2));
ft4 = reshape(resh(4,:),sizer(1)/4,sizer(2));

for ts = 1:nfreq*ncomp*ngraph
    graphresults{Rindstart+4,1} = 'mean';
    graphresults{Rindstart+5,1} = filetypes{1};
    graphresults{Rindstart+6,1} = filetypes{2};
    graphresults{Rindstart+7,1} = filetypes{3};
    graphresults{Rindstart+8,1} = filetypes{4};
    [h,p] = ttest2(resultsmat_mean(healthy,ts),resultsmat_mean(patients,ts),'Vartype','unequal');
    graphresults{Rindstart+4,ts+1} = p;
    [h,p] = ttest2(ft1(healthy,ts),ft1(patients,ts),'Vartype','unequal');
    graphresults{Rindstart+5,ts+1} = p;
    [h,p] = ttest2(ft2(healthy,ts),ft2(patients,ts),'Vartype','unequal');
    graphresults{Rindstart+6,ts+1} = p;
    [h,p] = ttest2(ft3(healthy,ts),ft3(patients,ts),'Vartype','unequal');
    graphresults{Rindstart+7,ts+1} = p;
    [h,p] = ttest2(ft4(healthy,ts),ft4(patients,ts),'Vartype','unequal');
    graphresults{Rindstart+8,ts+1} = p;
end

save('graphresults','graphresults');
