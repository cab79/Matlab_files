%make a violin plot using the violin FEX function
left = [1:2 5:6 9:10];
right = [3 4 7 8 11 12];
xlabel = {'','10% change','    ','30% change','    ','50% change'};
cmap = [
    0 1 0;
    1 1 0;
    0 0 1;
    0 1 1;
    1 0 0;
    1 0 1;
    ];

figure
violin(x(:,left),'xlabel', xlabel,'facecolor',cmap,'medc','k','mc','','bw',0.06)
ylabel('Response Time (s)','FontSize',14)
legend off
box off
