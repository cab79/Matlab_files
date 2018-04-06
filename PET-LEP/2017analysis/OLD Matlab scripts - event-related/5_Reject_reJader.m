clear all
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};


excludes = {
    '';
    '[2,31,33,34,35,36,39,41,46,47,51,52,58,59;]';
    
    
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

