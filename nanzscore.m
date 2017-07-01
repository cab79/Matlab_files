function zs = nanzscore(X,FLAG,DIM)

zs = bsxfun(@rdivide, bsxfun(@minus, X, nanmean(X,DIM)), nanstd(X,FLAG,DIM));