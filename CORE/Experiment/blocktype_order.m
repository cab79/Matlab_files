% blocktype order
close all
dname = 'C:\Users\cab79\Google Drive\2. Cambridge University\programs\Digitimer programs via labjack\Design templates';

fname = 'TSOT_design_template_part2_twohands.mat';

[hand,dc,cp,bi] = blocktype(dname,fname);

%% Counterbalance


%% plots

cmap = [
    [0.5 0.5 1] % L, 1DC, 10PC
    [1 0.25 0.25] % L, 3DC, 10PC
    [0.5 0.5 1] % R, 1DC, 10PC
    [1 0.25 0.25] % R, 3DC, 10PC
    [0.25 0.25 0.75] % L, 1DC, 30PC
    [0.75 0 0] % L, 3DC, 30PC
    [0.25 0.25 0.75] % R, 1DC, 30PC
    [0.75 0 0] % R, 3DC, 30PC
    [0 0 0.5] % L, 1DC, 50PC
    [0.5 0 0] % L, 3DC, 50PC
    [0 0 0.5] % R, 1DC, 50PC
    [0.5 0 0] % R, 3DC, 50PC
    [0.75 0.75 0.75]; %SideL
    [0.25 0.25 0.25]; %SideR
    ];

%figure
all = [hand;dc;cp;bi];
%imagesc(all)
%colormap(hsv)

% figure
% left = all(:,hand==1);
% [~,~,left(4,:)] = unique(left(4,:));
% imagesc(left)
% colormap(cmap)
% 
% figure
% right = all(:,hand==2);
% [~,~,right(4,:)] = unique(right(4,:));
% imagesc(right)
% colormap(cmap)

hv=unique(hand);
for hi=1:length(hv)
    all(1,all(1,:)==hv(hi))=max(all(4,:))+hi;
end
imagesc(all([1 4],:));
colormap(cmap)
set(gca, 'YTick', [],'FontSize',10);
xlabel('block')