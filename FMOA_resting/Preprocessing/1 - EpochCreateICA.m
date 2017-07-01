clear all
ds=2;
epoch_len = 4; %sec
% define subject numbers
subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2','OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

% define files to load
for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};

    s=[subject '_Base.cnt'];
    files=dir(s);
    


% load data from each file  
for f=1:length(files)  

    filename = files(f).name;
    tempcnt= loadcnt(filename);
    
%correct bad channels
if strcmp(num2str(subject), 'H8') == 1
    tempcnt.data(14,:) = mean(tempcnt.data([15 4 39 48 49 50 35 37],:),1);
    tempcnt.data(20,:) = mean(tempcnt.data([15 21 51 54 49 50 55 58],:),1);
    tempcnt.data(60,:) = mean(tempcnt.data([23 28],:),1);
elseif strcmp(num2str(subject), 'OA10') == 1
    tempcnt.data(27,:) = mean(tempcnt.data([21 58],:),1);
    tempcnt.data(61,:) = mean(tempcnt.data([27 29],:),1);
elseif strcmp(num2str(subject), 'OA12') == 1
    tempcnt.data(50,:) = mean(tempcnt.data([15 8 39 51],:),1);
elseif strcmp(num2str(subject), 'OA14') == 1
    tempcnt.data(18,:) = mean(tempcnt.data([11 24],:),1);
    tempcnt.data(63,:) = mean(tempcnt.data([21 62],:),1);
elseif strcmp(num2str(subject), 'OA16') == 1
    tempcnt.data(16,:) = mean(tempcnt.data([8 22 51 52],:),1);
    tempcnt.data(17,:) = mean(tempcnt.data([10 23],:),1);
elseif strcmp(num2str(subject), 'OA17') == 1
    tempcnt.data(16,:) = mean(tempcnt.data([8 22 51 52],:),1);
end

tempcnt.data = downsample(tempcnt.data',ds)'; % downsample to 125Hz
srate = tempcnt.header.rate/ds;

%hicutoff = 120;
lowcutoff = 0.5;
%filtorder = 6*fix(srate/hicutoff);
%tempcnt.data = eegfilt(tempcnt.data,srate,0,hicutoff);
tempcnt.data = eegfilt(tempcnt.data,srate,lowcutoff,0);

% extract triggers 
event = tempcnt.event;
stim = struct2cell(event);

% extract data point of triggers
offset = squeeze(stim(5,:,:));
offset=round([offset{:}]/ds);
stimtype = squeeze(stim(1,:,:));
stimtype=[stimtype{:}];

% find the index of each trigger value
ii1=find(stimtype==1);
ii2=find(stimtype==2);

if length(ii1) < 1
    error('no trigger 1 value')
end

if length(ii2) < 3
    error('not enough trigger 2 values')
end

% epoch, detrend and baseline correct

data1 = [];
data2 = [];

for i=1
    data1(:,:,1)=tempcnt.data(:, offset(ii1(i)):offset(ii1(i))+srate*60-1);
end

for i=2:3
    data2(:,:,i-1)=tempcnt.data(:, offset(ii2(i-1)):offset(ii2(i-1))+60*srate-1);
end

for i=4
    if size(tempcnt.data,2) < offset(ii2(i-1))+srate*60-1
        fileend = size(tempcnt.data,2);
    else
        fileend = offset(ii2(i-1))+srate*60-1;
    end
    endtemp = [tempcnt.data(:, offset(ii2(i-1)):fileend) zeros(size(tempcnt.data,1),offset(ii2(i-1))+srate*60-1-fileend)];
    data1(:,:,2)=endtemp(:,1:srate*60);
end

data1 = reshape(data1,size(data1,1),size(data1,2)*size(data1,3));
data2 = reshape(data2,size(data2,1),size(data2,2)*size(data2,3));

data_1 = [];
data_2 = [];

trial_number = (floor((size(data1,2))/srate)/epoch_len)-1;
for i=1:trial_number
    data_1(:,:,i)=data1(:, 1+epoch_len*srate*(i-1):epoch_len*srate*(i-1)+epoch_len*srate);
end
for i=1:size(data1,3)
    data_1(:,:,i)=(detrend2(squeeze(data_1(:,:,i))',100,100))';
end
data_1 = blcorrect4(data_1, 1);

trial_number = (floor((size(data2,2))/srate)/epoch_len)-1;
for i=1:trial_number
    data_2(:,:,i)=data2(:, 1+epoch_len*srate*(i-1):epoch_len*srate*(i-1)+epoch_len*srate);
end
for i=1:size(data2,3)
    data_2(:,:,i)=(detrend2(squeeze(data_2(:,:,i))',100,100))';
end
data_2 = blcorrect4(data_2, 1);


clear ii1 ii2 data1 data2

end % end of files loop

x=size(data_1);
eval(['save ' num2str(subject) '_data_matrix_dim1'  ' x']);
x=size(data_2);
eval(['save ' num2str(subject) '_data_matrix_dim2'  ' x']);

total_data = cat(3, data_1,data_2);
x=size(total_data);
Nelectrodes = x(1);
Nsamples = x(2);
Nevents = x(3);
total_data = reshape(total_data,Nelectrodes,Nevents*Nsamples);
eval(['save ' num2str(subject) '_data_epochs.mat'  ' total_data']);
b= jader(total_data, 30);
s=['b_' num2str(subject)];
save(s, 'b');


clear tempcnt total_data data_1 data_2 b

end % end of subjects loop