
% random number generaotr setting
rng(14,'twister');

% Clustering Fisher's Iris Data Using K-Means Clustering
load fisheriris
[cidx2,cmeans2] = kmeans(meas,3,'dist','sqeuclidean','display','iter'); % desired number of clusters set to 2, and using squared Euclidean distance
[silh2,h] = silhouette(meas,cidx2,'sqeuclidean'); % silhouette plot displays a measure of how close each point in one cluster is to points in the neighboring clusters

%the fourth measurement in these data, the petal width, is highly correlated with the third measurement and so a 3-D plot of the first three measurements gives a good 
%representation of the data, without resorting to four dimensions. COULD ALSO USE PCA COMPONENTS? If you plot the data, using different symbols for each cluster created by kmeans, you can identify 
%the points with small silhouette values, as those points that are close to points from other clusters.

ptsymb = {'bs','r^','md','go','c+'};
for i = 1:2
    clust = find(cidx2==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
% The centroids of each cluster are plotted using circled X's. 
plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko');
plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'kx'); 
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on

% unfinished - see http://uk.mathworks.com/help/stats/examples/cluster-analysis.html