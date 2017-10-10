function [Xsub,idx,ucol]=licols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
if ~nnz(X) %X has no non-zeros and hence no independent columns
     Xsub=[]; idx=[];
     return
end

if nargin<2, tol=1e-10; end
[Q, R, E] = qr(X,0); 
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);   
end

%Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
idx=sort(E(1:r));
Xsub=X(:,idx);    

% for the remainder, which indices are they most like?
[~,sorti]=sort(E);
ucol = [1:r nan(1,length(diagr)-r)];
idx_R = r+1:length(diagr);
for ir = idx_R 
   rval = corr(R(:,ir),R(:,1:r));
   [~,maxi] = max(rval);
   ucol(ir)=maxi;
end
ucol = ucol(sorti); % unique column indices of X