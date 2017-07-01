clear all

rawpth=('M:\Matlab\ExpStudy\Behaviour\Raw');
pth=('M:\Matlab\ExpStudy\Behaviour\Results');
condpth=('M:\Matlab\ExpStudy\Behaviour\Results\cond numbers');
cd(pth);

ncolhead=3;
[data hdr raw]=xlsread(fullfile(rawpth,'Pain_ratings_HGF.xls'));

nsub = length(hdr)-ncolhead;

intall = data(:,find(ismember(hdr,'int'))); 
uncertall = data(:,find(ismember(hdr,'uncert')));
condall = data(:,find(ismember(hdr,'cond')));
orderall = [1:60;61:120;121:180;181:240;241:300;301:360];

pint = [3 5 7 3 5 7];

for s=45%1:nsub
    %sort([1:16],'descend')
    
    % load inputs
    subname = hdr{ncolhead+s}
    prn=data(:,ncolhead+s);
    
    
    % Discrepent subjects only. After checking against output from
    % "TrialMatch.m"
    if strcmp(subname,'F14')
        block_order=[1 2 2 4:6]; % first block missing, second block is repetition of the first block
    else
        block_order=[1:6];
    end
    order_idx = reshape(orderall(block_order,:)',360,1);
    int = intall(order_idx);
    uncert = uncertall(order_idx);
    cond = condall(order_idx);
    
    %if isfield(results,subname)
    %    continue
    %end
    
    % find first and last non-NaN
    pr_idx0 = find(~isnan(prn));
    if ~isempty(pr_idx0); 
        pr_idx = pr_idx0([1,end]);
    else
        pr_idx = [1,length(prn)];
    end
    pr_idx1 = pr_idx(1):pr_idx(2);
    
    % data for analysis
    pr = prn(pr_idx1);
    int = int(pr_idx1);
    uncert = uncert(pr_idx1);
    cond = cond(pr_idx1);
    
    % remove trials from extended periods of no data
    prin = isnan(pr);
    prin_idx=zeros(length(prin),1);
    for i = 3:length(prin)-3
        if (prin(i) && prin(i+1) && prin(i+2)) || prin(i) && prin(i-1) && prin(i-2)
            prin_idx(i,1)=1;
        end
    end
    pr_idx1(find(prin_idx))=[];
    pr(find(prin_idx))=[];
    int(find(prin_idx))=[];
    uncert(find(prin_idx))=[];
    cond(find(prin_idx))=[];
    
    condsavename = fullfile(condpth,[subname '_conds.mat']);
    save(condsavename);
    
end