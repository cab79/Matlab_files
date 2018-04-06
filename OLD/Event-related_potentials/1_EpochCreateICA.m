clear all
Pdata ='I:\PET study\EEG session';
cd(Pdata)

ele = [1:30 33:64];
tr_start  = -2000;
tr_end = 749;
    
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};
%need: S4, P9, S15

% define files to load
for sub = length(subjects)-1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};

    s=[subject '*part2.vhdr'];
    eegfiles=dir(s);
    

% load data from each file  
%for f=1:length(eegfiles)  

filename = eegfiles.name;
tempcnt = pop_loadbv(Pdata, filename);
     
%correct bad channels
if strcmp(num2str(subject), 'P20_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P24_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P30_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P32_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P33_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P19_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
elseif strcmp(num2str(subject), 'P22_') == 1
    tempcnt.data(61,:) = mean(tempcnt.data([2 40],:),1);
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
elseif strcmp(num2str(subject), 'P23_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
elseif strcmp(num2str(subject), 'P27_') == 1
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
elseif strcmp(num2str(subject), 'P31_') == 1
    tempcnt.data(53,:) = mean(tempcnt.data([11 39],:),1);
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(33,:) = mean(tempcnt.data([17 3],:),1);
    tempcnt.data(34,:) = mean(tempcnt.data([17 4],:),1);
elseif strcmp(num2str(subject), 'P25_') == 1
    tempcnt.data(53,:) = mean(tempcnt.data([11 39],:),1);
    tempcnt.data(61,:) = mean(tempcnt.data([2 40],:),1);
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(33,:) = mean(tempcnt.data([17 3],:),1);
    tempcnt.data(34,:) = mean(tempcnt.data([17 4],:),1);
elseif strcmp(num2str(subject), 'P35_') == 1
    tempcnt.data(53,:) = mean(tempcnt.data([11 39],:),1);
    tempcnt.data(61,:) = mean(tempcnt.data([2 40],:),1);
    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
    tempcnt.data(33,:) = mean(tempcnt.data([17 3],:),1);
    tempcnt.data(34,:) = mean(tempcnt.data([17 4],:),1);
end

event = tempcnt.event;
stim = struct2cell(event);
chans = tempcnt.chanlocs;
chans= struct2cell(chans);
chans = squeeze(chans(1,:));
offset = squeeze(stim(1,:,:));
offset=[offset{:}];

stimtype = squeeze(stim(5,:,:))';
ii1=find(strcmp(stimtype,'S  1'));
ii2=find(strcmp(stimtype,'S  2'));

data1 = [];
data2 = [];

for i=1:length(ii1)
    data1(:,:,i)=(detrend(squeeze(tempcnt.data(:,offset(ii1(i))+tr_start:offset(ii1(i))+tr_end))'))';
end
data1 = blcorrect4(data1, 250);

for i=1:length(ii2)
    data2(:,:,i)=(detrend(squeeze(tempcnt.data(:,offset(ii2(i))+tr_start:offset(ii2(i))+tr_end))'))';
end
data2 = blcorrect4(data2, 250);

total_data = cat(3, data1,data2);

x=[size(total_data,1) size(total_data,2) length(ii1) length(ii2)];
eval(['save ' num2str(subject) 'data_matrix_dim'  ' x']);

if exist('total_data','var')
    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3)+x(4);
    total_data = reshape(total_data,Nelectrodes,Nevents*Nsamples);
    total_data = total_data(ele,:,:);
    eval(['save ' num2str(subject) 'data_segments.mat'  ' total_data']);
    b= jader(total_data, 30);
    s=[num2str(subject) 'b.mat'];
    save(s, 'b');
end


clear tempcnt total_data data1 data2

end % end of subjects loop