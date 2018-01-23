function W = gplotprepare_xlsxdata(S)
[NUM,TXT,RAW]=xlsread(fullfile(S.path,S.fname));
W.x = [RAW{S.nhead+1:end,find(strcmp(RAW(1,:),S.xvalhead))}]';
for p = 1:length(S.yvalhead)
    W(p).x = W(1).x;
    W(p).y = [RAW{S.nhead+1:end,find(strcmp(RAW(1,:),S.yvalhead{p}))}]';
    if isfield(S,'condhead')
        W(p).cond = [RAW{S.nhead+1:end,find(strcmp(RAW(1,:),S.condhead{p}))}]';
    end
    if isfield(S,'yvalhead')
        W(p).ptitle = S.yvalhead{p};
    end
end