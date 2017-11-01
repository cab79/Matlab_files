function [orth_wf,ucol,i] = find_ortho_symm_v3(S,wfall,tol_ord)

if ~S.ana_singlesub
    % concatenate over subjects and conditions
    wfall = {cat(2,wfall{:})};
end

for w = 1:length(wfall)
    wfall{w} = wfall{w}';

    % remove means over time
    wfall{w} = demean(wfall{w}, 1);
end


orth_wf = cell(size(wfall));
sizewf = size(wfall{1},1);
wfall = cat(1,wfall{:});

tryagain=1;
i = 0;
while tryagain
    i = i+1;
    % create rank-reduced from licols for concat
    [ld_wf,ld_sub,ucol,rank]=licols(wfall,tol_ord(i));
    disp(['tol = ' num2str(tol_ord(i)) '; rank = ' num2str(rank)])
    % de-concat
    ld_wf=reshape(ld_wf',size(ld_wf,2),sizewf,[]);
    %ld_wf = permute(ld_wf,[1 3 2]);
    ld_wf=squeeze(num2cell(ld_wf,[1 2]));
    % run each through ortho
    for w = 1:length(ld_wf)
        ld_wf{w} = ld_wf{w}';
        [orth_wf{w}, ~, ~, W,r,asize] = symmetric_orthogonalise_CAB(ld_wf{w}, 1);
        %orth_wf{w} =  orth_wf{w}';
        if isempty(orth_wf{w}) || r>S.outputrank
            % if any are rank deficient, start again with higher tol
            tryagain=1;
            %tol=tol*S.ranktolfactor;
            %disp(['try again on wf ' num2str(w) '/' num2str(length(ld_wf))])
            break
        elseif w == length(ld_wf)
            tryagain=0;
        end
    end
end
         
for w = 1:length(orth_wf)
    % remove means over time
    orth_wf{w} = demean(orth_wf{w}, 1);

    orth_wf{w} =  orth_wf{w}';
end
