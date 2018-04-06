
%exclude components
clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};
for x=1:length(subjects)

    subject = subjects(x);
    subject = char(subject);


    load ([subject 'data_segments.mat']);

    load ([subject 'b' ]);


    excludes = {
       '[1,3,9,10,11,12,15,17,20;]'
       '[1,2,3,4,6,7,8,10,28;]'
       '[1,2,7,12,16,27;]'
       '[1,2,3,4,6;]'
       '[1,2,3,4,6,8,9;]'
       '[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,21;]'
       '[1,2,4,5,8,9,19;]'
       '[1,3,7,8;]'
       '[1,11;]'
       '[1,13;]'
       '[1,2,3,4,16;]'
       '[1,2,6;]'
       '[1,2,3,16;]'
       '[1,2,5,7,8,9,12,13,15;]'
       '[1,4,6;]'
       '[1,2,5;]'
       '[1,2,4,5,6,7,13,14,19;]'
       '[1,2,4,8,15;]'
       '[1,2,4,5,6,11,12;]'
       '[1,3,4,6,8,10,11,13,14,15,16,17,20,21,22,23,24,25,26,27,29,30;]'
       '[1,2,3,5,6,8,9,10,11,14,21,23,24,25,26,28,30;]'
       '1 2 3 5 7 8 11 24'
       '1 2 3 5 6 7 8 19 21 22'
       '1 2 3 11 12 21 27'
       '1 2 3 8'
       '1 3 5 7 8'
       '1 2 3 4 5 6 10 25 26 27'
       '1 2 3 5 7 9 10'
       '1 2 3 4 10 18'
       '1 2 5 6 7'
       '1 2 14 21'
       '1 4 5'
       '1 2 3 4 7 8 9 30'
    };


    exclude=excludes(x);
    exclude=[exclude{:}];

    exclude=str2num(exclude);

    include = [1:30];
    include(exclude) = [];


    act = b * total_data;
    ib = pinv(b);

    iact= ib(:, include) * act(include,:);  

    %iact = iirfilt(iact,500,0,30); %filter the data

    Nelectrodes = 62;
    Nsamples = 2750;
    Nevents = size(iact,2)/Nsamples;

    total_data_ICA=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

    eval(['save ' subject 'total_data_ICA.mat'  ' total_data_ICA']);
end

clear exclude include total_data total_data_ICA

