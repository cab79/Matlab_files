clear all
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};


excludes = {
    '13 16 28 39 43 46 49'
    '18 24 16 35 39 41 56'
    '2 17 57 58'
    '3 4 10 11 14 15 31 33 36 41 42 56 57'
    '5 6 10 15 23 27 31 37 40 41 54'
    '6:14 50:53 60' % P8 check ICA for slow waves
    '11 13 16 19 35 37 45:47'
    '' % P16 check ICA for HR
    '3 4 15 34 50 58'
    '32'
    '2 10 12 21 29 39 42 44 59'
    '12 19 21 41' %S5 check ICA for blinks at end of trial
    '' %S6 check ICA for slow waves
    '19 23 28 43 59' 
    '8 12 13 17 18 23 24 32 34 49 57' %S9 check ICA
    '25 40 45 59' %S10 check ICA for slow waves
    '9 19 22 31 45 46'
    '31'
    '1:3 21 31 41'
    '1 2 6 15 16 19 20 26 27 30 31 32 54 55 58'
    '8 14 15 17 20 22 27 28 35 38 41 51'
    '10 18 25 33 46'
    '59 60'
    '1 2 8 9 10 14:17 18 19 23 24 28 30:33 40 42 43 44 47 56'
    '1 2 4 5 6 7 11 12 19 23 28 32 46 56 57'
    '16 31 46'
    '1 16 47 48'
    '1 7 16 29 39 40'
    '2 3 15 23 36 45 37 49:51 53'
    '7 8 17 38 '
    '5 55 59'
    '1 2 17 21 26 31 46'
    '6 15 16 18 26 28 35 53 55 58:60'
};

for sub = length(subjects):length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    filename = 'total_data_ICA';
    load ([subject filename '.mat']);
    
    load([subject 'data_matrix_dim']);

    reject = excludes(sub);
    reject=reject{:};
    reject=str2num(reject);
    
    tt1 = 1:x(3);
    tt2 = x(3)+1:x(3)+x(4);
    ex1 = find(ismember(tt1,reject));
    ex2 = find(ismember(tt2,reject));
    n1 = x(3)-length(ex1);
    n2 = x(4)-length(ex2);
 
    total_data_ICA(:,:,reject) = [];
    
    x = size(total_data_ICA);
    x(3:4) = [n1 n2];
    eval(['save ' num2str(subject) 'data_matrix_dim2 x']);

    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3)+x(4);
    total_data_ICA = reshape(total_data_ICA,Nelectrodes,Nevents*Nsamples);
    eval(['save ' num2str(subject) 'data_segments2.mat total_data_ICA']);
    b= jader(total_data_ICA, 30);
    s=[num2str(subject) 'b2.mat'];
    save(s, 'b');

    clear b s total_data_ICA x reject tt1 tt2 n1 n2

end

