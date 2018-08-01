function [D] = HGF_traj2mat(D,S,opt)
% D is a struct with length subjects containing an HGF field with
% trajectories
% S contains S.traj 

for i = 1:length(D)
    D(i).trajmat=[];
    D(i).trajmat_group=[];
    D(i).trajmat_label={};
    % for each traj group... reformat into double array
    for tg = 1:length(S.traj)
        for pm = 1:size(S.traj{tg})
            pm_name = S.traj{tg}{pm,1};
            for tr = 1:length(S.traj{tg}{pm,2})
                tr_name = S.traj{tg}{pm,2}(tr);
                tr_dat = D(i).HGF.fit.traj.(pm_name{:}).(tr_name{:});
                D(i).trajmat = [D(i).trajmat, tr_dat];
                nvar=size(tr_dat,2);
                D(i).trajmat_group = [D(i).trajmat_group, tg*ones(1,nvar)];
                for nv = 1:nvar
                    D(i).trajmat_label = [D(i).trajmat_label, {[pm_name{:} '_' tr_name{:} num2str(nv)]}];
                end
            end
        end
    end
    
    if isfield(opt,'PCA') && opt.PCA==1
       
        % determine numcomponent by doing an eig on the covariance matrix
        covar = D(i).trajmat'*D(i).trajmat;
        [PC, diag] = eig(covar);
        [Diag,ind] = sort(diag(Diag),'descend');
        Diag = Diag ./ sum(Diag);
        Dcum = cumsum(Diag);
        numcomp = find(Dcum>.9999,1,'first');
        D(i).PCdiag = Diag(1:numcomp);
        D(i).PC = PC(1:numcomp);
    end
    
end
    
