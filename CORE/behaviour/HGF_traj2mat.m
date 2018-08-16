function [S] = HGF_traj2mat(S,D)
% D is a struct with length subjects containing an HGF field with
% trajectories
% S contains S.traj 

% format HGF data
% for each traj group... reformat into predictor matri
S.pred=[];
S.pred_group=[];
S.pred_label={};
for tg = 1:length(S.traj)
    for pm = 1:size(S.traj{tg})
        pm_name = S.traj{tg}{pm,1};
        for tr = 1:length(S.traj{tg}{pm,2})
            tr_name = S.traj{tg}{pm,2}(tr);
            tr_dat = D.HGF.fit.traj.(pm_name{:}).(tr_name{:});
            if ~isempty(S.traj{tg}{pm,4}{tr})
                % select levels
                tr_dat = tr_dat(:,S.traj{tg}{pm,4}{tr});
            end
            if S.traj{tg}{pm,3}{1}(tr)
                tr_dat = abs(tr_dat);
            end
            S.pred = [S.pred, tr_dat];
            nvar=size(tr_dat,2);
            S.pred_group = [S.pred_group, tg*ones(1,nvar)];
            for nv = 1:nvar
                try
                    lev = S.traj{tg}{pm,4}{tr}(nv);
                catch
                    lev=nv;
                end
                S.pred_label = [S.pred_label, {[pm_name{:} '_' tr_name{:} num2str(lev)]}];
            end
        end
    end
end
