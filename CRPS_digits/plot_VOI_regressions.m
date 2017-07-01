%% 
X = [];%X axis regressor
Y = [];%VOI

%%
close all
nGrps = [13 13];
nConds = 10;
nPnts = sum(nGrps)*nConds;

r =[];
pr =[];
t =[];
pt =[];
C = colormap(hsv(max(nGrps)));
for g = 1:length(nGrps)
    figure(g)
    for s = 1:nGrps(g)
        start_ind = (g-1)*130 + (s-1)*nConds + 1;
        end_ind = (g-1)*130 + (s)*nConds;
        ind = start_ind:end_ind;
        x=zscore(X(ind));
        y=zscore(Y(ind));
        scatter(x,y,'MarkerFaceColor', C(s,:));
        lsline
        [r(g,s) pr(g,s)] = corr(x,y,'type','Pearson');
        stats = regstats(y,x);
        t(g,s) = stats.tstat.t(2);
        pt(g,s) = stats.tstat.t(1);
        hold on
    end
end

figure(g+1)
scatter(ones(1,length(r')),r(1,:))
hold on
scatter(2*ones(1,length(r')),r(2,:))

figure(g+2)
scatter(ones(1,length(t')),t(1,:))
hold on
scatter(2*ones(1,length(t')),t(2,:))

