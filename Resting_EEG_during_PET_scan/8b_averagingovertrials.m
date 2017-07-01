clear all

subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_'};
sessions = {'2','3'};
scans = {'pain','nonpain'};

for sess = 1:length(scans)
    session = sessions(sess);
    session = session{:};

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
load([char(subjects(sub)) 'total_data_ICA' session '_acc.mat']);
load([subject 'data_dim_acc' session]);
eval(['data = total_data_ICA' session ';']);

if x(4)+x(5)+x(6) ~= size(data,3); errordlg('trials per condition do not match data size');end
t1 = 1:x(4);
t2 = x(4)+1:x(4)+x(5);
t3 = x(4)+x(5)+1:x(4)+x(5)+x(6);
data1 = data(:,:,t1);
data2 = data(:,:,t2);
data3 = data(:,:,t3);

clear data t1 t2 t3 x

for dt = 1:3
    
eval(['data = data' num2str(dt) ';' ]);

for i=1:size(data,3)
    data(:,:,i)=(detrend2(squeeze(data(:,:,i))',250,250))';
end

data = blcorrect4(data, 1);

avg = squeeze(mean(data,3));

name=[subject 'avg_' scans{sess} '_' num2str(dt)];
save(name, 'avg'); clear avg;

end
end
end


