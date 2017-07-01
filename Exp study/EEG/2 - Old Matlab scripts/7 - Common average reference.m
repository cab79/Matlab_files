subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

for n = 1:length(subjects)

fnames={'_1_avgP2.mat';
    '_2_avgP2.mat';
    '_3_avgP2.mat';
    '_4_avgP2.mat';
    '_5_avgP2.mat';
    '_6_avgP2.mat';
  };

for x=1:length(fnames);

fname=char(fnames(x));
subject=char(subjects(n));
fname= [subject fname(1:8)];
load(fname);
mcr= mean(avg([1:2 4:30 33:64],:));
for i=1:64
    avg(i,:)=avg(i,:)- mcr;
end
eval(['save ' fname '_ca.mat'  ' avg'])

end

fnames={'_1_avg.mat';
    '_2_avg.mat';
    '_3_avg.mat';
    '_4_avg.mat';
    '_5_avg.mat';
    '_6_avg.mat';
  };

for x=1:length(fnames);

fname=char(fnames(x));
subject=char(subjects(n));
fname= [subject fname(1:6)];
load(fname);
mcr= mean(avg([1:2 4:30 33:64],:));
for i=1:64
    avg(i,:)=avg(i,:)- mcr;
end
eval(['save ' fname '_ca.mat'  ' avg'])

end

end