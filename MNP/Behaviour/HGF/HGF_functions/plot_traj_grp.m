
level = 2;
sub_ind = [ones(16,1)' 2*ones(15,1)' 3*ones(16,1)'];
col = {'k','r','b'};
figure
for g = 1:3
    hold all
    trajgrp = squeeze(mu_traj(find(sub_ind==g),level,:));
    %tmedian = nanmedian(trajgrp,1);
    tmean = nanmean(trajgrp,1);
    lower = tmean-nanstd(trajgrp,1);
    upper = tmean+nanstd(trajgrp,1);

    ts = 1:length(tmean);
    
    plot(upper, col{g}, 'LineWidth', 1);
    plot(lower, col{g}, 'LineWidth', 1);
    fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
         col{g}, 'EdgeAlpha', 0, 'FaceAlpha', 0.15);

end