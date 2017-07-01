clear all

% INPUT Args:
subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};
% reduced n
%subjects = {'M1';'M2';'M6';'M7';'M10';'M14';'M15';'M16';'M17';'M29';'M32';'M33';'M20';'M21';'M22';'M24';'M25';'M26';'M35';'M36';'M37';'M38';'M40'};

NSub = length(subjects);

elect_extract = 'FCZ.'; % OR use:
%chan_extract = [46 18 45 12 13 14 51 25 19]; %9 ele around Cz

conditions = 2;
results = zeros(NSub,6);

%%%%%%%%%%
fnames={'_1_p2_av.dat'
    '_2_p2_av.dat'
    '_3_p2_av.dat'
    '_4_p2_av.dat'
    '_5_p2_av.dat'
    '_6_p2_av.dat'
  };

if exist('elect_extract') ==1
[ele,x,y,chans] = textread('C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan.locs2','%f %f %f %s');
chans = upper(chans);
chan_extract = find(strcmp(chans, elect_extract));
end

for n = 1:NSub
sub = subjects(n);
c1 = load([char(sub) char(fnames(1))]);
c2 = load([char(sub) char(fnames(2))]);
c3 = load([char(sub) char(fnames(3))]);
c4 = load([char(sub) char(fnames(4))]);
c5 = load([char(sub) char(fnames(5))]);
c6 = load([char(sub) char(fnames(6))]);

c1 = mean(c1(1,chan_extract),2);
c2 = mean(c2(1,chan_extract),2);
c3 = mean(c3(1,chan_extract),2);
c4 = mean(c4(1,chan_extract),2);
c5 = mean(c5(1,chan_extract),2);
c6 = mean(c6(1,chan_extract),2);

results(n,1) = c1;
results(n,2) = c2;
results(n,3) = c3;
results(n,4) = c4;
results(n,5) = c5;
results(n,6) = c6;

end

% normalise the data to percent of mean across conditions, for each subject
for i=1:NSub
    mean_results= mean(results(i,1:4),2);
    %results(i,1:4) = ((results(i,1:4)-mean_results)/mean_results)*100;
    results(i,1:4) = results(i,1:4)-mean_results;
    mean_results= mean(results(i,5:6),2);
    %results(i,5:6) = ((results(i,5:6)-mean_results)/mean_results)*100;
    results(i,5:6) = results(i,5:6)-mean_results;
end
