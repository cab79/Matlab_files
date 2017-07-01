clear all
if isunix
    filepath = '/scratch/cb802/Data/CIP_resting/fMRI';
else
    filepath = 'W:\Data\CIP_resting\fMRI';
end
cd(filepath);
conntype = 'fMRI_pearsons';
anatype = '';
cd(fullfile(filepath,conntype,anatype));
files = dir(fullfile(filepath,conntype,anatype,'*_graph.mat'));
healthy = 1:39;
%healthy = [1 2 3 4 5 20 24 32 34]; % age and sex matched
patients = 40;
allsub = 1:length([healthy patients]);
%allsub = 1:length(files);
files = files(allsub);
filetypes = {''};
nfreq = 1;
ncomp = 3;
ngraph = 8;
%nodROI = [70 76]; % L and R post insula
nodROI = []; % loop through all nodes if blank
%nodROI = [26 28]; % L and R ant insula
%nodROI = [85 87 94 105 112 132]; % L and R precuneus
nonROI = 1:160;
tval = 2; % t value 
load age;

clear graphresults
graphresults{1,1} = 'clustering';
graphresults{1,2} = 'characteristic path length';
graphresults{1,3} = 'global efficiency';
graphresults{1,4} = 'modularity';
graphresults{1,5} = 'modules';
graphresults{1,6} = 'centrality';
graphresults{1,7} = 'participation coefficient';
graphresults{1,8} = 'connection density';

graphresults = reshape(repmat(graphresults',1,nfreq*ncomp)',1,nfreq*ncomp*ngraph);
subhead1{1,1} = '';
%subhead1{1,2} = 'theta';
%subhead1{1,3} = 'alpha';
%subhead1{1,4} = 'beta';
%subhead1{1,5} = 'gamma';
subhead1 = repmat(reshape(repmat(subhead1',1,ncomp)',1,nfreq*ncomp),1,ngraph);
graphresults(2,:) = subhead1;

subhead2{1,1} = 'allnodes';
subhead2{1,2} = 'ROI';
subhead2{1,3} = 'ROI_ratio';
subhead2 = repmat(subhead2,1,nfreq*ngraph);
graphresults(3,:) = subhead2;

graphresults(:,2:(nfreq*ncomp*ngraph+1)) = graphresults(:,1:(nfreq*ncomp*ngraph));
graphresults(1:3,1) = cell(3,1);

if isempty(nodROI)
    nd = nonROI;
else 
    nd = nodROI;
end

out_range = nan(160,size(graphresults,2)-1,length(files));
all_range = out_range;

for n = nd
    n
    nROI = nonROI;
    nROI(n) = [];
    for f = 1:length(files)
        for ft = 1:length(filetypes)
            Rindstart = (f-1)*length(filetypes) + (ft-1) + 1;
            filename = files(f).name;
            [pth basename ext] = fileparts(filename);
            graphresults{3+Rindstart,1} = [basename(1:10) filetypes{ft}];
            datfile = dir(fullfile(filepath,conntype,anatype,[basename(1:10) '*' filetypes{ft} '*_graph.mat']));
            randfile = dir(fullfile(filepath,conntype,anatype,[basename(1:10) '*' filetypes{ft} '*_randgraph.mat']));
            dat = load(datfile.name);
            rand = load(randfile.name);
            for gt = 1:ngraph
                clear datamat
                Cindstart = nfreq*ncomp*(gt-1)+1;
                gdat = dat.graphdata{gt,2};
                rdat = rand.graphdata{gt,2};
                %if ~isreal(gdat) 
                %    gdat = abs(gdat);
                %    rdat = abs(rdat);
                %end
                rdat = mean(rdat,length(size(rdat))); % 
                rdat = squeeze(mean(rdat,1));
                gdat = squeeze(mean(gdat,1));
                normdat = gdat./rdat; 
                %normdat = normdat(isfinite(normdat));
                %fin = isfinite(normdat);
                %in_ind = find(fin==1);
                if size(normdat,2)>1
                    datamat(1,:) = nanmean(normdat(:,nROI),2)';
                    datamat(2,:) = nanmean(normdat(:,n),2)';
                    datamat(3,:) = datamat(2,:)./datamat(1,:);
                else 
                    datamat = [normdat';nan(1,nfreq);nan(1,nfreq)];
                end
                datamat = reshape(datamat,1,nfreq*ncomp);
                resultsmat(Rindstart,Cindstart:Cindstart+(ncomp*nfreq-1)) = datamat;
            end
        end
    end
    
    
    for col = 1:size(resultsmat,2)
        meancol = mean(resultsmat(:,col));
        [B,BINT,R] = regress(resultsmat(:,col),age(allsub));
        %scatter(resultsmat(:,col),meancol+R)
        resultsmat(:,col) = meancol+R;
    end
    
    graphresults(4:3+Rindstart,2:nfreq*ncomp*ngraph+1) = num2cell(resultsmat);
    
    % for each subject/file, calculate how no. of SDs of their data outside
    % the mean, and whether their values of out of range.
    for f = 1:length(files)
        cas = f;
        con = allsub;
        con(f) = [];
        for col = 1:size(resultsmat,2)
            se = std(resultsmat(con,col))/sqrt(size(resultsmat(con,col),1));
            sepi = std(resultsmat(con,col))*sqrt(1 + 1/size(resultsmat(con,col),1));
            av = mean(resultsmat(con,col));
            
            CIthresh = av+tval*se
            PIthresh = av+tval*sepi
            
            if ((resultsmat(f,col)>max(resultsmat(con,col)) && resultsmat(f,col) > PIthresh) || (resultsmat(f,col)<min(resultsmat(con,col)) && resultsmat(f,col) < -PIthresh))
                out_range(n,col,f) = (resultsmat(f,col) - av)/sepi;
            end
            if (resultsmat(f,col) > PIthresh || resultsmat(f,col) < -PIthresh)
                all_range(n,col,f) = (resultsmat(f,col) - av)/sepi;
            end
        end
    end
end

save('PI_tval_data_ageregressed','out_range','all_range','graphresults');

%for ts = 1:nfreq*ncomp*ngraph
%    %[h,p] = ttest2(resultsmat(healthy,ts),resultsmat(patients,ts),'Vartype','unequal');
%    %graphresults{Rindstart+4,ts+1} = p;
%    graphresults{Rindstart+4,ts+1} = mean(resultsmat(healthy,ts));
%end

%save('graphresults','graphresults');

%compROI = [26 28 44 48 55 56 59 61 70 76]; % ant, post and mid insula
compROI = [14 19 21 27];
%compROI = [76]; % left post insula
nout = nan(size(all_range,3),size(all_range,2)); % number of regions*parameters that are out of range for that file/subject and parameter
for c = 1:size(all_range,2) % parameter
    vplot = [];
    for f = 1:size(all_range,3) % file/subject
        vmat = [];
        outnn = ~isnan(squeeze(out_range(compROI,c,f)));
        if any(outnn)
            nout(f,c) = length(find(outnn==1));
            vout = squeeze(out_range(compROI(find(outnn==1)),c,f));
            vmat(:,2) = vout;
            vmat(:,1) = f*ones(length(vout),1);
            vplot = [vplot; vmat];
        end
    end
    if size(vplot,1)>0
        figure;
        scatter(vplot(:,1),vplot(:,2));
        title([graphresults{1,1+c} '__' graphresults{3,1+c}]);
    end
end

close all
for c = 1:size(all_range,2)
    vplot = [];
    for f = 1:size(all_range,3)
        vmat = [];
        outnn = ~isnan(squeeze(all_range(compROI,c,f)));
        if any(outnn)
            %nout(f,c) = length(find(outnn==1));
            vout = squeeze(all_range(compROI(find(outnn==1)),c,f));
            vmat(:,2) = vout;
            vmat(:,1) = f*ones(length(vout),1);
            vplot = [vplot; vmat];
        end
    end
    if size(vplot,1)>0
        figure;
        scatter(vplot(:,1),vplot(:,2));
        title([graphresults{1,1+c} '__' graphresults{3,1+c}]);
    end
end

patient_out_range = out_range(:,:,40);

