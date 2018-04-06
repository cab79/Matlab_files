clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P25_';'P31_'};

for dt = 1:2

    grand_avg=zeros(62,2750);

    for sub = 1:length(subjects)
        subject = subjects(sub);
        subject = subject{:};

        load([subject 'avg_' num2str(dt) '_ca.mat'])
        grand_avg=grand_avg + avg;
        clear avg       
    end
    grand_avg = grand_avg / sub; 
    eval(['save Grand_avg_' num2str(dt) '_ca.mat'  ' grand_avg']);
    clear grand_avg
end





