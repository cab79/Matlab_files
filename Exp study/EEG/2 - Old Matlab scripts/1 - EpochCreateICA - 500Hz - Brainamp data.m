clear all

Pdata = 'C:\Documents and Settings\mdmoscab\Desktop\Chris Data\Expectancy Study';
subjects = {'OA16'}; 
ele = [1:2 4:30 33:64];
%dsample = 10;
tr_start  = 2000;
tr_end = 749;

for sub = 1
    subject = subjects(sub);
    subject = subject{:};

    cd(Pdata)
    filename=[subject '.vhdr'];
    
    tempcnt = pop_loadbv(Pdata, filename);
    %tempcnt.data = downsample(tempcnt.data',dsample)';

% extract triggers
%stim = tempcnt.event.stimtype;

event = tempcnt.event;
stim = struct2cell(event);

chans = tempcnt.chanlocs;
chans= struct2cell(chans);
chans = squeeze(chans(1,:));

%correct reference and bad channels
tempcnt.data(65,:) = mean(tempcnt.data([17 18 21 22],:),1);
chans(1,65) = {'FCz'};
tempcnt.data(66,:) = mean(tempcnt.data([45 59],:),1);
chans(1,66) = {'PO5'};
tempcnt.data(67,:) = mean(tempcnt.data([46 60],:),1);
chans(1,67) = {'PO6'};

if strcmp(num2str(subject), 'EF6_1') == 1
    tempcnt.data(33,:) = mean(tempcnt.data([34 40 42 48],:),1); %F4
%    tempcnt.data(30,:) = mean(tempcnt.data([64 62],:),1);
end

example_cnt= loadcnt('OA15.cnt');
example_chans = cell(1,length(example_cnt.electloc));
for c = 1:length(example_cnt.electloc)
    example_chans(1,c) = {example_cnt.electloc(1,c).lab};
end
chans = upper(chans);
example_chans = upper(example_chans);
chans(find(strcmp('TP9',chans) == 1)) = {'M1'};
chans(find(strcmp('EOG',chans) == 1)) = {'HEOG'};
chans(find(strcmp('ECG',chans) == 1)) = {'VEOG'};

tempcnt2 = tempcnt;
tempcnt2.data = tempcnt2.data(1:64,:);
for c = 1:length(example_chans)
    e = find(strcmp(example_chans(c),chans) == 1);
    tempcnt2.data(c,:) = tempcnt.data(e,:);
end
tempcnt = tempcnt2;

offset = squeeze(stim(1,:,:));
offset=[offset{:}];
offset = offset/10;
stimtype = squeeze(stim(5,:,:))';
S127s = strcmp('S127',stimtype);
stimtype(S127s) = [];
boundary = strcmp('boundary',stimtype);
stimtype(boundary) = [];

if strcmp(num2str(subject), 'EF6_1') == 1
    load([subject '_offsets500Hz.mat']);
    load stimtype
    stimtype = stimtype(2:113);
    %tr_start  = 4000;
    %tr_end = 2749;
end

ii1=find(strcmp(stimtype,'S 63'));
ii2=find(strcmp(stimtype,'S 95'));
ii3=find(strcmp(stimtype,'S 31'));
ii4=find(strcmp(stimtype,'S111'));
ii5=find(strcmp(stimtype,'S 47'));
ii6=find(strcmp(stimtype,'S 79'));

if exist('ii1') == 1
trial_number = length(ii1);
for i=1:trial_number
    if offset(ii1(i))+749 < size(tempcnt.data,2)
    data1(:,:,i)=tempcnt.data(:, offset(ii1(i))-tr_start:offset(ii1(i))+tr_end);
    end
end
if ii1>0
for i=1:size(data1,3)
    data1(:,:,i)=(detrend(squeeze(data1(:,:,i))'))';
end
data1 = blcorrect4(data1, 250);
end
end

if exist('ii2') == 1
    trial_number = length(ii2);
for i=1:trial_number
    if offset(ii2(i))+tr_end < size(tempcnt.data,2) && offset(ii2(i))-tr_start > 0
    data2(:,:,i)=tempcnt.data(:, offset(ii2(i))-tr_start:offset(ii2(i))+tr_end);
    end
end
if ii2>0
for i=1:size(data2,3)
    data2(:,:,i)=(detrend(squeeze(data2(:,:,i))'))';
end
data2 = blcorrect4(data2, 250);
end
end

if exist('ii3') == 1 
    trial_number = length(ii3);
for i=1:trial_number
    if offset(ii3(i))+tr_end < size(tempcnt.data,2)
    data3(:,:,i)=tempcnt.data(:, offset(ii3(i))-tr_start:offset(ii3(i))+tr_end);
    end
end
if ii3>0
for i=1:size(data3,3)
    data3(:,:,i)=(detrend(squeeze(data3(:,:,i))'))';
end
data3 = blcorrect4(data3, 250);
end
end

if exist('ii4') == 1 
    trial_number = length(ii4);
for i=1:trial_number
    if offset(ii4(i))+tr_end < size(tempcnt.data,2)
    data4(:,:,i)=tempcnt.data(:, offset(ii4(i))-tr_start:offset(ii4(i))+tr_end);
    end
end
if ii4>0
for i=1:size(data4,3)
    data4(:,:,i)=(detrend(squeeze(data4(:,:,i))'))';
end
data4 = blcorrect4(data4, 250);
end
end

if exist('ii5') == 1 
    trial_number = length(ii5);
for i=1:trial_number
    if offset(ii5(i))+tr_end < size(tempcnt.data,2) && offset(ii5(i))-tr_start > 0
    data5(:,:,i)=tempcnt.data(:, offset(ii5(i))-tr_start:offset(ii5(i))+tr_end);
    end
end
if ii5>0
for i=1:size(data5,3)
    data5(:,:,i)=(detrend(squeeze(data5(:,:,i))'))';
end
data5 = blcorrect4(data5, 250);
end
end

if exist('ii6') == 1 
    trial_number = length(ii6);
for i=1:trial_number
    if offset(ii6(i))+tr_end < size(tempcnt.data,2) && offset(ii6(i))-tr_start > 0
    data6(:,:,i)=tempcnt.data(:, offset(ii6(i))-tr_start:offset(ii6(i))+tr_end);
    end
end
if ii6>0
for i=1:size(data6,3)
    data6(:,:,i)=(detrend(squeeze(data6(:,:,i))'))';
end
data6 = blcorrect4(data6, 250);
end
end

total_data = cat(3, data1, data2, data3, data4, data5, data6);
total_data = total_data(ele,:,:);

stim_matrix = [length(ii1) length(ii2) length(ii3) length(ii4) length(ii5) length(ii6)];
eval(['save ' num2str(subject) '_stim_matrix'  ' stim_matrix']);

x=size(total_data);
eval(['save ' num2str(subject) '_data_matrix_dim'  ' x']);


Nelectrodes = x(1);
Nsamples = x(2);
Nevents = x(3);

total_data = reshape(total_data,Nelectrodes,Nevents*Nsamples);

eval(['save ' num2str(subject) '_data_epochs.mat'  ' total_data']);

b= jader(total_data, 30);
s=['b_' num2str(subject)];
save(s, 'b');





clear ii1 ii2 ii3 ii4 ii5 ii6 total_data data1 data2 data3 data4 data5 data6

end % end of subjects loop









