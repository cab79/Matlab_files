for t=500:10:749 %%time from : to every 10 data points
t2=10+t
topoplot(mean(avg([1:2 4:30 33:64], t:t2),2), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'), title(t); %%change if chan.locs2 is in a different filw
   pause;	
end

for t=1:30 %%time from : to every 10 data points
topoplot(b(t, [1:2 4:30 33:64]), 'C:\Documents and Settings\All Users\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2'), title(t); %%change if chan.locs2 is in a different filw
   pause;	
end