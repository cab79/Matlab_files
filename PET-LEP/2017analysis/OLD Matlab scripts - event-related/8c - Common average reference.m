
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};


for n = length(subjects):length(subjects)
    fnames={'avg_1.mat';
        'avg_2.mat';
      };
    for x=1:length(fnames);
        fname=char(fnames(x));
        subject=char(subjects(n));
        fname= [subject fname(1:5)];
        load(fname);
        mcr= mean(avg([1:62],:));
        for i=1:62
            avg(i,:)=avg(i,:)- mcr;
        end
        eval(['save ' fname '_ca.mat'  ' avg'])
    end
end