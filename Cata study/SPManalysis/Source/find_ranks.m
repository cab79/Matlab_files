function [tol_ord,S] = find_ranks(S,wfall,inputclus,clusdir)

if ~S.ana_singlesub
    % concatenate over subjects and conditions
    wfall = {cat(2,wfall{:})};
end

for w = 1:length(wfall)
    wfall{w} = wfall{w}';

    % remove means over time
    wfall{w} = demean(wfall{w}, 1);
end

% METHOD:
% concatenate
sizewf = size(wfall{1},1);
wfall = cat(1,wfall{:});
tol=5e-9; % MAY NEED TO INCREASE TO PREVENT POSITIVE DEFINITE ERROR LATER. Probably caused by low-pass filtering, which reduces the rank again.

rank=3; % any number greater than 2 to get the loop going
S.rankinfo = [];
i = 0;
while rank>2
    i = i+1;
    % create rank-reduced from licols for concat
    [~,~,ucol,rank]=licols(wfall,tol);
    
    % create new cluster img and find number of connected regions within each
    % region
    img = create_reduced_image(S,inputclus,clusdir,0,ucol);
    [ucom,IA,IB] = unique(img);
    nReg=nan(1,length(ucom));
    for u = 1:length(ucom)
        if ucom(u)==0
            continue
        end
        % how many connected areas?
        binimg = zeros(size(img));
        binimg(IB==u) = 1;
        CC = bwconncomp(binimg,6);
        nReg(u)=length(CC.PixelIdxList);
    end
    nReg_av = nanmean(nReg);
    
    % calcuate score to minimise
    S.rankinfo(i,1) = rank;
    S.rankinfo(i,2) = nReg_av;
    if strcmp(S.rankminmax,'min')
        S.rankinfo(i,3) = rank*nReg_av^2;
    elseif strcmp(S.rankminmax,'max')
        S.rankinfo(i,3) = (1/rank)*nReg_av^2;
    end
    S.rankinfo(i,4) = tol;
    
    % display
    disp(['i = ' num2str(i) '; tol = ' num2str(tol) '; rank = ' num2str(rank) '; nReg_av = ' num2str(nReg_av)])
    
    % set up for next iteration
    tol=tol*S.ranktolfactor;
end

% create list of tols in order of preference
[~,ord] = sort(S.rankinfo(:,3));
S.rankinfo = S.rankinfo(ord,:);
tol_ord = S.rankinfo(:,4);
