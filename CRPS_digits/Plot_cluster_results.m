clear results_mat
load clustersERP
load Acc_cov_leftright

for s = 1:length(results)
    for s2 = 1:length(results{1,s}) 
        for r = 1:size(results{1,s}{1,s2},2)-1
            results_mat(r,s,s2) = results{1,s}{1,s2}{2,r+1};
        end
    end
end
results_mat = reshape(permute(results_mat,[1 3 2]),r,s*s2);
%results_12 = results_mat1-results_mat2;

close all
ind = {1:13;14:26;27:39;40:52};
R=[];P=[];
for i = 1:4
    figure(i);
    x=cov(ind{i});
    y=results_mat(1,ind{i});
    scatter(x,y)
    [r p] = corrcoef(x,y);
    R(i)=r(1,2);
    P(i)=p(1,2);
end
R
P