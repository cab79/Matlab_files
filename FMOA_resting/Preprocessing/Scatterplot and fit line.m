clear all
row_start = 3; % row data starts
data_lines = 1:15; % which rows of data to use
plim = 0.1;

[data] = textread('Behavioural - ERP correlations 3.txt','','delimiter','   ','headerlines',1);
data = data(data_lines,row_start:size(data,2));
no_fact = size(data,2);

nan_pos = isnan(data);
nmean = nanmean(data);

for f = 1:size(data,2)
    r = find(nan_pos(:,f) == 1);
    data(r,f) = nmean(f);
end

[R,P]=corrcoef(data(:,:));

fact = textread('Behavioural - ERP correlations 3.txt','%s');
fact = fact(row_start:no_fact+row_start-1);
fact_no = 1:length(fact);

for f = 2:no_fact
    sig = find(P(1:f-1,f) < plim);
    results(f,1) = fact(f);
    results(f,2:1+length(sig)) = fact(sig);
end

results_fact = fact(sig)
results_pval = results(fact_no(sig))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fact1 = 1:no_fact;
fact2 = 1:no_fact;

% go through all correlations
f=0;
for b = 44
for e = 65:1000
    e = e+f;
A = data(:,fact2(e));
B = data(:,fact1(b));
corr = R(fact2(e),fact1(b));
pval = P(fact2(e),fact1(b));
name = [fact(fact2(e)) '_____' fact(fact1(b)) '_____' num2str(corr) '_____' num2str(pval)];
plot(A,B,'bo'), title(name);
lsline
ch = getkey;
if ch == 108
    continue
elseif ch == 107
    f = f-2;
end
end
end


% go through significant correlations only
for f = 65:no_fact; %factor number (row of 'results')
sig = find(P(1:f-1,f) < plim);
for c = 2:length(sig) % number of sig correlate (column of 'results')
A = data(:,f);
B = data(:,sig(c-1));
corr = R(sig(c-1),f);
pval = P(sig(c-1),f);
name = [results{f,1} '_____' results{f,c} '_____' num2str(corr) '_____' num2str(pval)];
plot(A,B,'bo'), title(name);
lsline
pause
end
end

% non-parametric test
f = 87;
c = 3;
%d = 2;
A = data(:,f);
B(:,1) = data(:,sig(c-1));
%B(:,2) = data(:,sig(c-1));
[betanonp tstats Fnp] = nonparamReg(A,B);
