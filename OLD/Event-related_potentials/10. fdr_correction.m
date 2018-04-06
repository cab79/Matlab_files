% load data
filename = 'p_values.csv';
p = load(filename);

% FDR correction
q=0.05;
[pID,pN] = fdr(p,q);

if any(pN)
    msgbox('significant results found');
    sig_results = find(p>=pID)
else
    msgbox('no significant results');
end