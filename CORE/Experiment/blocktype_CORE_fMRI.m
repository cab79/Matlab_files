function [mode,cp,bi,blockii,btypes] = blocktype_CORE_fMRI(dname,fname)

load(fullfile(dname,fname));

btypes = {
    [1 2],1,1 % Tactile, 20Pc
    [3 4],2,1 % Auditory, 20Pc
    [5 6],1,2 % Tactile, 40Pc
    [7 8],2,2 % Auditory, 40Pc
};
cond = seq.condnum;
cond(1)=0; % first trial is a new block
blockii = find(cond==0);
blocks = cond(1,blockii+1);

mode=[];
cp=[];
bi=[];
for b = 1:length(blocks)
    for bt = 1:length(btypes(:,1))
        if ismember(blocks(b),btypes{bt,1})
            bi(b)=bt;
            break
        end
    end
    mode(b)=[btypes{bi(b),2}];
    cp(b)=[btypes{bi(b),3}];
end
