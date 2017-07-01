subject = 'OA2';
load([subject '_total_data_ica.mat']);

load ([subject '_data_matrix_dim']);
Nelectrodes = x(1);
Nsamples = x(2);
Nevents = x(3);

total_data_ICA = reshape(total_data_ICA,Nelectrodes,Nevents*Nsamples);
b= jader(total_data_ICA, 30);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

act = b * total_data_ICA; 
ib = pinv(b);
eegplot(act); 