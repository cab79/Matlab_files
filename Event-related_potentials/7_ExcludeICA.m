
%exclude components
clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};

for i=length(subjects):length(subjects)

    subject = subjects(i);
    subject = char(subject);


    load ([subject 'data_segments2.mat']);

    load ([subject 'b2' ]);


    excludes = {
       '[11;]'
       ''
       ''
       ''
       '[4,30;]'
       '[1,4,9,30;]'
       '[30;]'
       ''
       ''
       '25'
       ''
       '6'
       '[1,2,6,17;]'
       '5'
       '[2,8;]'
       '[3,6,14;]'
       '9'
       ''
       ''
       ''
       '[2,17,30;]'
       ''
       ''
       ''
       ''
       ''
       ''
       ''
       '2'
       '6'
       '5'
       '2'
       '1 2 3 8 9'
    };


    exclude=excludes(i);
    exclude=[exclude{:}];

    exclude=str2num(exclude);

    include = [1:30];
    include(exclude) = [];


    act = b * total_data_ICA;
    ib = pinv(b);

    iact= ib(:, include) * act(include,:);  

    %iact = iirfilt(iact,500,0,30); %filter the data
    
    %for da = 1:size(iact,2)
    %    mcr= mean(iact(:,da),1);
    %    iact(:,da) = iact(:,da) - mcr;
    %end

    Nelectrodes = 62;
    Nsamples = 2750;
    Nevents = size(iact,2)/Nsamples;

    total_data_ICA2=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

    eval(['save ' subject 'total_data_ICA2.mat'  ' total_data_ICA2']);
end

clear exclude include total_data_ICA total_data_ICA2

