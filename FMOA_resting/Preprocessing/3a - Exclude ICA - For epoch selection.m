clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2','OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

for s=1:length(subjects)

subject = subjects(s);
subject = char(subject);

    load ([subject '_data_epochs.mat']);
    load (['b_' subject]);


excludes = {
'[1,2,3,8;]'%1
'[1,2,4;]'%2
'[1,2;]'%3
'[1,2,3,4,17;]'%4
'[1,2,3;]'%5
'[1,16;]'%6
'[1,2;]'%7
'[1,2,14;]'%8
'[1,2,4,5,6;]'%9
'[1,2,3,6;]'%10
'[1,2,4,5,6,7;]'%11
'[1;]'%12
'[1,2,3;]'%13
'[1,2,3,4,6;]'%14
'[1,4,6,7;]'%15
'[1,2,5;]'%16
'[1,2,4;]'
'[1,2,3,5,7;]'
'[1,2,9,11,18;]'
'[1,2,3,5;]'%4
'[1,2,5,7,11;]'%5
'[4;]'%6
'[1,2;]'%7
'[1,4;]'%8
'[1,2,8;]'%9
'[1,3,17;]'%10
'[1,9,18;]'%11
'[2;]'%12
'[1,7,14;]'%13
'[1,4,6,18;]'%14
'[1,2,3,9;]'%15
'[1;]'%16
'[1,2,11,15;]'%1
'[1,2,4,6;]'%2
'[1;]'%4
'[1,2,3,4,12;]'%5
'[1,15;]'%6
'[1,2;]'%7
'[1,7,25;]'%8
'[1,2,3,4,5;]'%9
'[1,4;]'%10
'[1,5;]'%11
'[1,2,6;]'%12
'[1,3;]'%13
'[1,6,7;]'%14
'[1,8,9,10;]'%15
'[1,5;]'%16
'[1,3,5;]'%17
};


exclude=excludes(s);
exclude=[exclude{:}];

exclude=str2num(exclude);

include = [1:30];
include(exclude) = [];


act = b * total_data;
ib = pinv(b);

iact= ib(:, include) * act(include,:);  

load ([subject '_data_matrix_dim1']);
d1 = x;
load ([subject '_data_matrix_dim2']);
d2 = x;
clear x
x = [d1(1) d1(2) d1(3)+d2(3)];

Nelectrodes = x(1);
Nsamples = x(2);
Nevents = x(3);

total_data_ICA=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

eval(['save ' subject '_total_data_ICA.mat'  ' total_data_ICA']);


clear exclude include total_data total_data_ICA x

end