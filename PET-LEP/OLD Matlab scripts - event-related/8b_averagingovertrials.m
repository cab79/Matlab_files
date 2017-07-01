clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
    load([char(subjects(sub)) 'total_data_ICA2.mat']);
    load([subject 'data_matrix_dim2.mat']);
    eval(['data = total_data_ICA2;']);

if x(3)+x(4) ~= size(data,3); errordlg('trials per condition do not match data size');end
t1 = 1:x(3);
t2 = x(3)+1:x(3)+x(4);
data1 = data(:,:,t1);
data2 = data(:,:,t2);

clear data t1 t2 x

for dt = 1:2
    
eval(['data = data' num2str(dt) ';' ]);

for i=1:size(data,3)
    data(:,:,i)=(detrend(squeeze(data(:,:,i))'))';
end

data = blcorrect4(data, 250);

avg = squeeze(mean(data,3));

name=[subject 'avg_' num2str(dt)];
save(name, 'avg'); clear avg;

end
end



