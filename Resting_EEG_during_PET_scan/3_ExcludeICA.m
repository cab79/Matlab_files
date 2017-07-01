
%exclude components
clear all

subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
for x=1:length(subjects)

subject = subjects(x);
subject = char(subject);



    load ([subject 'data_segments2.mat']);

    load ([subject 'b2' ]);


excludes = {
   '' %S1 
   ''
   ''
   ''
   '' %S9 big impedence?
   ''
   ''
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

for da = 1:size(iact,2)
    mcr= mean(iact(:,da),1);
    iact(:,da) = iact(:,da) - mcr;
end

Nelectrodes = 62;
Nsamples = 2500;
Nevents = size(iact,2)/Nsamples;

total_data_ICA2=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

eval(['save ' subject 'total_data_ICA2.mat'  ' total_data_ICA2']);
end

clear exclude include total_data total_data_ICA2


clear all
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
for x=1:length(subjects)
subject = subjects(x);
subject = char(subject);

    load ([subject 'data_segments3.mat']);

    load ([subject 'b3' ]);

%excludes = {
%   '2'
%   '3 4'
%   '1 2 6 12'
%   '1 4 5'
%   '2 3 4'
%   '1 3'
%};
excludes = {
   '' %S1
   ''
   ''
   ''
   ''
   ''
   ''
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

for da = 1:size(iact,2)
    mcr= mean(iact(:,da),1);
    iact(:,da) = iact(:,da) - mcr;
end

Nelectrodes = 62;
Nsamples = 2500;
Nevents = size(iact,2)/Nsamples;

total_data_ICA3=reshape(iact, Nelectrodes, Nsamples, Nevents); % re-froms the epochs

eval(['save ' subject 'total_data_ICA3.mat'  ' total_data_ICA3']);


clear exclude include total_data total_data_ICA3 x

end