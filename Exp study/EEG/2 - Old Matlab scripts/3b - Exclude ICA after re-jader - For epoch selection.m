clear all;

% type in subject numbers

subjects = {'H3', 'F1', 'F4', 'F7', 'OA2'};

for i=1:length(subjects)

subject = subjects(i);
subject = char(subject);



    load ([subject '_total_data_ica.mat']);

    load (['b' subject]);
    
    load ([subject '_data_matrix_dim']);

    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3);
    
    total_data_ICA = reshape(total_data_ICA,Nelectrodes,Nevents*Nsamples);


excludes = {'4 11'
    '11 17'
    '4'
    '1 2 3'
    '2 4 5'
};


exclude=excludes(i);
exclude=[exclude{:}];

exclude=str2num(exclude);

include = [1:30];
include(exclude) = [];


act = b * total_data_ICA;
ib = pinv(b);

iact= ib(:, include) * act(include,:);  




total_data_ICA2=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

eval(['save ' subject '_total_data_ica.mat'  ' total_data_ICA2']);


clear exclude include total_data_ICA2 total_data_ICA x

end