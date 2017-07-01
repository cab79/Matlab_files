%Check Averages
clear all
scans = {'pain','nonpain'};
t=1:1501;

col = {'r','b'};

for dt = 1:3

    for sess = 1:length(scans)
    
        load(['Grand_avg_' scans{sess} '_' num2str(dt) '.mat'])
        plot(t, -grand_avg(18,:), col{sess})
        hold on

    end
    hold off
    pause
end