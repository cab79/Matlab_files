function [hand,dc,cp,bi,blockii,btypes] = blocktype(dname,fname)

load(fullfile(dname,fname));

btypes = {
    [1 3],1,1,1 % L, 1DC, 10PC
    [2 4],1,2,1 % L, 3DC, 10PC
    [5 7],2,1,1 % R, 1DC, 10PC
    [6 8],2,2,1 % R, 3DC, 10PC
    [9 11],1,1,2 % L, 1DC, 30PC
    [10 12],1,2,2 % L, 3DC, 30PC
    [13 15],2,1,2 % R, 1DC, 30PC
    [14 16],2,2,2 % R, 3DC, 30PC
    [17 19],1,1,3 % L, 1DC, 50PC
    [18 20],1,2,3 % L, 3DC, 50PC
    [21 23],2,1,3 % R, 1DC, 50PC
    [22 24],2,2,3 % R, 3DC, 50PC
};
if exist('dt','var'); design=dt.design; end;
design(2,1)=0; % first trial is a new block
blockii = find(design(2,:)==0);
blocks = design(2,blockii+1);

hand=[];
dc=[];
cp=[];
bi=[];
for b = 1:length(blocks)
    for bt = 1:length(btypes(:,1))
        if ismember(blocks(b),btypes{bt,1})
            bi(b)=bt;
            break
        end
    end
    hand(b)=[btypes{bi(b),2}];
    dc(b)=[btypes{bi(b),3}];
    cp(b)=[btypes{bi(b),4}];
end
