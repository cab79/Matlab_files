clear all
close all
dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\export_fig')
addpath('C:\Data\CORE\eeg')
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

%% figure options
fontsize = 12; 
save_figs=1;
xticks = [-200:4:800];
chan_time = 'extent_time';
clusdir='Grp_clusters';
plotclus = {'c1_spm','c2_spm','c3_spm','c4_spm','c5_spm','c6_spm','c7_spm','c8_spm','c9_spm','c10_spm','c11_spm','c12_spm','c13_spm','c14_spm','c15_spm'};
subplotgap=0.01;
xaxisgap = 0.6;
yaxisgap=0.15;

%% SPM clusters


i=1;
S.con(i).spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_stats_BRR_all_chan_condHGF_notrans_20190221T154622_pred4_spm';
S.con(i).clusdir=clusdir;
S.con(i).plotclus = plotclus;
S.con(i).colormap = cbrewer('seq', 'Greens', 100, 'pchip');
S.con(i).dsample=1;
S.con(i).ylabel='epsi 3';

i=2;
S.con(i).spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_stats_BRR_all_chan_condHGF_notrans_20190221T154622_pred3_spm';
S.con(i).clusdir=clusdir;
S.con(i).plotclus = plotclus;
S.con(i).colormap = cbrewer('seq', 'Blues', 100, 'pchip');
S.con(i).dsample=1;
S.con(i).ylabel='epsi 2';

i=3;
S.con(i).spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_stats_BRR_all_chan_condHGF_notrans_20190221T154622_pred2_spm';
S.con(i).clusdir=clusdir;
S.con(i).plotclus = plotclus;
S.con(i).colormap = cbrewer('seq', 'Purples', 100, 'pchip');
S.con(i).dsample=1;
S.con(i).ylabel='epsi 1';

i=4;
S.con(i).spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_stats_BRR_all_chan_condHGF_notrans_20190221T154622_pred1_spm';
S.con(i).clusdir=clusdir;
S.con(i).plotclus = plotclus;
S.con(i).colormap = cbrewer('seq', 'YlOrBr', 100, 'pchip');
S.con(i).dsample=4;
S.con(i).ylabel='digit change';

i=5;
S.con(i).spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\sensor\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_spm';
S.con(i).clusdir=clusdir;
S.con(i).plotclus = plotclus;
S.con(i).colormap = cbrewer('seq', 'Greys', 100, 'pchip');
S.con(i).dsample=4;
S.con(i).ylabel='original';


figure
for i = 1:length(S.con)
    ax(i) = subplot('Position', [yaxisgap, 1-(((1-xaxisgap)*i)/length(S.con))+subplotgap, 1-yaxisgap, (1-xaxisgap)/length(S.con)-subplotgap]);
    hold on
    for c = 1:length(S.con(i).plotclus)
        try
            Cnii = load_nii(fullfile(S.con(i).spm_path,S.con(i).clusdir,[S.con(i).plotclus{c} '.nii']));
        catch
            continue
        end
        Cdim = size(Cnii.img);
        Cmat=reshape(permute(Cnii.img,[2 1 3]),Cdim(1)*Cdim(2),[]);
        switch chan_time
            case 'chan_time'
                
            case 'std_time'
                Cmat=nanstd(Cmat,[],1);
                
            case 'extent_time'
                Cmat=sum(Cmat>0);
        end
        tindex = 1:S.con(i).dsample:Cdim(3)-S.con(i).dsample;
        imagesc(xticks,[],Cmat(:,tindex),'AlphaData',Cmat(:,tindex)>0); 
    end
    colormap(ax(i),S.con(i).colormap);
    colorbar('Box','off')
    ylabel(S.con(i).ylabel,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
    cl(:,i) = caxis;
    if i~=length(S.con)
        set(ax(i), 'XColor',[1 1 1]);
    end
    set(ax(i),'YColor',S.con(i).colormap(90,:));
    hold off
end
set(ax, 'CLim', [mean(cl(1,:)), mean(cl(2,(cl(2,:)>1)))]);
set(ax, 'YTick', [],'Xlim',[min(xticks) max(xticks)],'Ylim',[0.48 1.52]);