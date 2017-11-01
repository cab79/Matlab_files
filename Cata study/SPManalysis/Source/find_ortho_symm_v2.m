function [orth_wf,ucol] = find_ortho_symm_v2(S,wfall,inputclus,clusdir)

orth_wf = cell(size(wfall));
% METHOD:
% concatenate
sizewf = size(wfall{1},1);
wfall = cat(1,wfall{:});
tol=5e-9; % MAY NEED TO INCREASE TO PREVENT POSITIVE DEFINITE ERROR LATER. Probably caused by low-pass filtering, which reduces the rank again.
tryagain=1;
while tryagain
    % create rank-reduced from licols for concat
    [ld_wf,ld_sub,ucol,rank]=licols(wfall,tol);
    disp(['tol = ' num2str(tol) '; rank = ' num2str(rank)])
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
            tol=tol*S.ranktolfactor;
            disp(['try again on wf ' num2str(w) '/' num2str(length(ld_wf))])
            break
        elseif w == length(ld_wf)
            %new_clus_img = create_reduced_image(S,inputclus,clusdir,0)
            %score = nReg*av_nClus +1;
            tryagain=0;
        end
    end
end
