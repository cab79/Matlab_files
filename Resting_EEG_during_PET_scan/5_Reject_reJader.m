clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};

excludes = {
    ''
    ''
    ''
    ''
    ''
    ''
    ''
};

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    filename = 'total_data_ICA2';
    load ([subject filename '.mat']);
    
    load([subject 'data_matrix_dim_2']);

    reject = excludes(sub);
    reject=reject{:};
    reject=str2num(reject);
    
    tt1 = 1:(x(3)/3);
    tt2 = (x(3)/3)+1:2*(x(3)/3);
    tt3 = 2*(x(3)/3)+1:3*(x(3)/3);
    ex1 = find(ismember(tt1,reject));
    ex2 = find(ismember(tt2,reject));
    ex3 = find(ismember(tt3,reject));
    n1 = (x(3)/3)-length(ex1);
    n2 = (x(3)/3)-length(ex2);
    n3 = (x(3)/3)-length(ex3);
 
    total_data_ICA2(:,:,reject) = [];
    
    x = size(total_data_ICA2);
    x(4:6) = [n1 n2 n3];
    eval(['save ' num2str(subject) 'data_dim_acc2 x']);

    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3);
    total_data_ICA2 = reshape(total_data_ICA2,Nelectrodes,Nevents*Nsamples);
    eval(['save ' num2str(subject) 'data_segments2_acc.mat total_data_ICA2']);
    b= jader(total_data_ICA2, 30);
    s=[num2str(subject) 'b2_acc.mat'];
    save(s, 'b');

    clear b s total_data_ICA2 x reject tt1 tt2 tt3

end

clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
excludes = {
    ''
    ''
    ''
    ''
    ''
    ''
    ''
};
 
 
for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    filename = 'total_data_ICA3';
    load ([subject filename '.mat']);
    
    load([subject 'data_matrix_dim_3']);

    reject = excludes(sub);
    reject=reject{:};
    reject=str2num(reject);
    
    tt1 = 1:(x(3)/3);
    tt2 = (x(3)/3)+1:2*(x(3)/3);
    tt3 = 2*(x(3)/3)+1:3*(x(3)/3);
    ex1 = find(ismember(tt1,reject));
    ex2 = find(ismember(tt2,reject));
    ex3 = find(ismember(tt3,reject));
    n1 = (x(3)/3)-length(ex1);
    n2 = (x(3)/3)-length(ex2);
    n3 = (x(3)/3)-length(ex3);
 
    total_data_ICA3(:,:,reject) = [];
    
    x = size(total_data_ICA3);
    x(4:6) = [n1 n2 n3];
    eval(['save ' num2str(subject) 'data_dim_acc3 x']);

    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3);
    total_data_ICA3 = reshape(total_data_ICA3,Nelectrodes,Nevents*Nsamples);
    eval(['save ' num2str(subject) 'data_segments3_acc.mat total_data_ICA3']);
    b= jader(total_data_ICA3, 30);
    s=[num2str(subject) 'b3_acc.mat'];
    save(s, 'b');

    clear b s total_data_ICA3 x reject tt1 tt2 tt3

    clear events

end
