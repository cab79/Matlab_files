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
            try
                tr_dat = D.HGF.fit.traj.(pm_name{:}).(tr_name{:});
            catch
                % derived traj
                if strcmp(tr_name{:},'surp')
                    m1h = D.HGF.fit.traj.(pm_name{:}).muhat(:,1);
                    u = D.HGF.u(1,:);
                    poo = m1h.^u.*(1-m1h).^(1-u); % probability of observed outcome
                    tr_dat = -log2(poo);
                elseif strcmp(tr_name{:},'inferv')
                    mu2 = D.HGF.fit.traj.(pm_name{:}).mu(:,2);
                    sa2 = D.HGF.fit.traj.(pm_name{:}).sa(:,2);
                    sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                    tr_dat = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
                elseif strcmp(tr_name{:},'pv')
                    mu2 = D.HGF.fit.traj.(pm_name{:}).mu(:,2);
                    mu3 = D.HGF.fit.traj.(pm_name{:}).mu(:,3);
                    sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                    tr_dat = sigmoid_mu2.*(1 -sigmoid_mu2).*exp(mu3); 
                end
                
            end
            if ~isempty(S.traj{tg}{pm,4}{tr})
                % select levels
                tr_dat = tr_dat(:,S.traj{tg}{pm,4}{tr});
            end
            if S.traj{tg}{pm,3}{1}(tr)==1 % absolute
                tr_dat = abs(tr_dat);
            elseif S.traj{tg}{pm,3}{1}(tr)==2 % rectified
                tr_dat(tr_dat<0) = 0;
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
