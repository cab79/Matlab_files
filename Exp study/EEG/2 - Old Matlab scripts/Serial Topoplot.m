load M2_2_avg_ca

figure
for t=500:688 %%time from : to every 10 data points
t2=2+t
topoplot(mean(grand_avg([1:2 4:30 33:64], t:t2),2), 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'), title(t); %%change if chan.locs2 is in a different filw
   pause;	
end
