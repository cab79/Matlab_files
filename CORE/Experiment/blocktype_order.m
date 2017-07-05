% blocktype order

dname = 'C:\Users\cab79\Google Drive\2. Cambridge University\programs\Digitimer programs via labjack\Design templates';

fname = 'TSOT_design_template_part4_twohands_reversed.mat';

[hand,dc,cp,bi] = blocktype(dname,fname);

%% Counterbalance


%% plots

cmap = [
    1 1 0;
    1 0 1;
    0 1 1;
    1 0 0;
    0 1 0;
    0 0 1;
    ];

%figure
all = [hand;dc;cp;bi];
%imagesc(all)
%colormap(hsv)

figure
left = all(:,hand==1);
[~,~,left(4,:)] = unique(left(4,:));
imagesc(left)
colormap(cmap)

figure
right = all(:,hand==2);
[~,~,right(4,:)] = unique(right(4,:));
imagesc(right)
colormap(cmap)