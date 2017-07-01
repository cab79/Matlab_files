function calcgraph(basename,conntype,mattype,varargin)

if isunix
    filepath = '/scratch/cb802/Data/CRPS_resting/EEG/';
else
    filepath = 'W:\Data\CRPS_resting\EEG\';
end
%loadsubj

param = finputcheck(varargin, {
    'randomise', 'string', {'on','off','edgerand'}, 'off'; ...
    'latticise', 'string', {'on','off'}, 'off'; ...
    'numrand', 'integer', [], 50; ...
    'rewire', 'integer', [], 50; ...
    'heuristic', 'integer', [], 50; ...
    });

load chandist
chandist = chandist / max(chandist(:)); % needed for calculation of modular span only.

tvals = 0.5:-0.025:0.1;
%tvals = 0.2:-0.01:0.1;

if strcmp(param.randomise,'on') || strcmp(param.randomise,'edgerand')
    %savename = sprintf('%s/%s/%s%srandgraph.mat',filepath,conntype,basename,conntype);
    savename = fullfile(filepath,conntype,[basename '_' conntype '_randgraph.mat']);
    numruns = param.numrand;
elseif strcmp(param.latticise,'on')
    distdiag = repmat(1:length(sortedlocs),[length(sortedlocs) 1]);
    for d = 1:size(distdiag,1)
        distdiag(d,:) = abs(distdiag(d,:) - d);
    end
    distdiag = distdiag ./ max(distdiag(:));
    %savename = sprintf('%s/%s/%s%slattgraph.mat',filepath,conntype,basename,conntype);
    savename = fullfile(filepath,conntype,[basename '_' conntype '_latgraph.mat']);
    numruns = param.numrand;
else
    %savename = sprintf('%s/%s/%s%sgraph.mat',filepath,conntype,basename,conntype);
    savename = fullfile(filepath,conntype,[basename '_' conntype '_graph.mat']);
    numruns = 1;
end

graphdata{1,1} = 'clustering';
graphdata{2,1} = 'characteristic path length';
graphdata{3,1} = 'global efficiency';
graphdata{4,1} = 'modularity';
graphdata{5,1} = 'modules';
graphdata{6,1} = 'centrality';
graphdata{7,1} = 'modular span';
graphdata{8,1} = 'participation coefficient';
graphdata{9,1} = 'connection density';
graphdata{10,1} = 'mutual information';

fprintf('Processing %s',basename);

load([filepath conntype filesep basename '_' conntype '.mat']); % contains matrix and bootmat
if ~isempty(mattype)
    eval(['matrix = matrix_' mattype ';']);
    eval(['bootmat = bootmat_' mattype ';']);
end

if strcmp(param.randomise,'edgerand')
    bootmat=[];
end

%[sortedchan,sortidx] = sort({chanlocs.labels});
%if ~strcmp(chanlist,cell2mat(sortedchan))
%    error('Channel names do not match!');
%end
%matrix = matrix(:,sortidx,sortidx);
%bootmat = bootmat(:,sortidx,sortidx,:);
%     pval = pval(:,sortidx,sortidx);

%chanlocs = chanlocs(sortidx);
%     chanXYZ = [cell2mat({chanlocs.X})' cell2mat({chanlocs.Y})' cell2mat({chanlocs.Z})'];

for f = 1:size(matrix,1)
    %if f == 2; keyboard;end
    for iter = 1:numruns % for each randomised matrix
        fprintf('ITER %d',iter);

        for thresh = 1:length(tvals)
            fprintf(' %d',thresh);
            if strcmp(param.randomise,'edgerand')
                %randomisation
                if ~isempty(bootmat)
                    cohmat = squeeze(bootmat(f,:,:,iter));
                else 
                    cohmat=[];
                end
            else
                cohmat = squeeze(matrix(f,:,:));
            end
            
            cohmat(isnan(cohmat)) = 0;
            cohmat = abs(cohmat);
            
            if isempty(bootmat) && isempty(cohmat)
                cohmat = squeeze(matrix(f,:,:));
                cohmat(isnan(cohmat)) = 0;
                cohmat = abs(cohmat);
                weicoh = threshold_proportional(cohmat,tvals(thresh));
                [weicoh eff] = randmio_und_connected(weicoh, 10);
                
                %eff
            else
                weicoh = threshold_proportional(cohmat,tvals(thresh));
            end
            %                 bincoh = double(threshold_proportional(cohmat,tvals(thresh)) ~= 0);
            
            %%%%%%  WEIGHTED %%%%%%%%%
            
            allcc{iter}(thresh,:) = clustering_coef_wu(weicoh);
            allcp{iter}(thresh) = charpath(distance_wei(weight_conversion(weicoh,'lengths')));
            alleff{iter}(thresh) = efficiency_wei(weicoh);
            allbet{iter}(thresh,:) = betweenness_wei(weight_conversion(weicoh,'lengths'));
            allden{iter}(thresh) = density_und(weicoh);
            
            for i = 1:param.heuristic
                [Ci, allQ{iter}(thresh,i)] = modularity_louvain_und(weicoh);
                
                allCi{iter}(thresh,i,:) = Ci;
                
                modspan = zeros(1,max(Ci));
                for m = 1:max(Ci)
                    if sum(Ci == m) > 1
                        distmat = chandist(Ci == m,Ci == m) .* weicoh(Ci == m,Ci == m);
                        distmat = nonzeros(triu(distmat,1));
                        modspan(m) = sum(distmat)/sum(Ci == m);
                    end
                end
                %iter = iter
                %thresh = thresh
                %i = i
                
                allms{iter}(thresh,i) = max(nonzeros(modspan));
                
                allpc{iter}(thresh,i,:) = participation_coef(weicoh,Ci);
            end
        end
    end
    
    for iter = 1:numruns
        for thresh = 1:length(tvals)
            %clustering coeffcient
            graphdata{1,2}(f,thresh,1:length(chanlocs),iter) = allcc{iter}(thresh,:);
            
            %characteristic path length
            graphdata{2,2}(f,thresh,iter) = allcp{iter}(thresh);
            
            %global efficiency
            graphdata{3,2}(f,thresh,iter) = alleff{iter}(thresh);
            
            % modularity
            graphdata{4,2}(f,thresh,iter) = mean(allQ{iter}(thresh,:));
            
            % community structure
            graphdata{5,2}(f,thresh,1:length(chanlocs),iter) = squeeze(allCi{iter}(thresh,1,:));
            
            %betweenness centrality
            graphdata{6,2}(f,thresh,1:length(chanlocs),iter) = allbet{iter}(thresh,:);
            
            %modular span
            graphdata{7,2}(f,thresh,iter) = mean(allms{iter}(thresh,:));
            
            %participation coefficient
            graphdata{8,2}(f,thresh,1:length(chanlocs),iter) = mean(squeeze(allpc{iter}(thresh,:,:)));
            
            %connection density
            graphdata{9,2}(f,thresh,iter) = allden{iter}(thresh);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %                 %BINARY
            %
            %                 for i = 1:param.heuristic
            %                     [Ci, allQ(i)] = modularity_louvain_und(bincoh);
            %
            %                     modspan = zeros(1,max(Ci));
            %                     for m = 1:max(Ci)
            %                         if sum(Ci == m) > 1
            %                             distmat = chandist(Ci == m,Ci == m) .* bincoh(Ci == m,Ci == m);
            %                             distmat = nonzeros(triu(distmat,1));
            %                             modspan(m) = sum(distmat)/sum(Ci == m);
            %                         end
            %                     end
            %                     allms(i) = max(nonzeros(modspan));
            %
            %                     allpc(i,:) = participation_coef(bincoh,Ci);
            %                 end
            %
            %                 %clustering coefficient
            %                 graphdata{1,3}(f,thresh,1:length(chanlocs),iter) = clustering_coef_bu(bincoh);
            %
            %                 %characteristic path length
            %                 graphdata{2,3}(f,thresh,iter) = charpath(distance_bin(bincoh));
            %
            %                 %global efficiency
            %                 graphdata{3,3}(f,thresh,iter) = efficiency_bin(bincoh);
            %
            %                 %modularity
            %                 graphdata{4,3}(f,thresh,iter) = mean(allQ);
            %
            %                 %community structure
            %                 graphdata{5,3}(f,thresh,1:length(chanlocs),iter) = Ci;
            %
            %                 %betweenness centrality
            %                 graphdata{6,3}(f,thresh,1:length(chanlocs),iter) = betweenness_bin(bincoh);
            %
            %                 %modular span
            %                 graphdata{7,3}(f,thresh,iter) = mean(allms);
            %
            %                 %participation coefficient
            %                 graphdata{8,3}(f,thresh,1:length(chanlocs),iter) = mean(allpc);
            %
            %                 %connection density
            %                 graphdata{9,3}(f,thresh,iter) = density_und(bincoh);
            %
            %                 %             %rentian scaling
            %                 %             [N, E] = rentian_scaling(bincoh,chanXYZ,5000);
            %                 %             E = E(N<size(bincoh,1)/2);
            %                 %             N = N(N<size(bincoh,1)/2);
            %                 %             b = robustfit(log10(N),log10(E));
            %                 %             graphdata{9,3}(f,thresh) = b(2);
            %
            %
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
end
fprintf('\n');

save(savename, 'graphdata', 'tvals');