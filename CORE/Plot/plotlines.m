fname = 'C:\Data\CORE\SPMstats\t-200_299_b-200_0_m_0_299_CP_Grp_Odd_Subject_4_cleaned_tm_spm\Odd_clusters\VOI_c1_short_20170804T185519.xlsx'
dat = xlsread(fname);
cols = [1 2];
ncols=3;
x = [1;2];
close all

% plot all conditions 
if 0
    for n=1:ncols
        y = dat(:,cols+2*(n-1))';
        figure
        line(x,y)
    end
end

% average and then plot
y=[];
for n = 1:length(cols)
    y(n,:) = mean(dat(:,cols(n)+[0 2 4]),2)';
end
figure
line(x,y)