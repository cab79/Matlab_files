clear all
Pdata ='C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG from PET sessions';
Presults = 'C:\Documents and Settings\mdmoscab\Desktop\Chris Data\PET DPN study\EEG from PET sessions\whole_scan';

ele = [1:30 33:64];
tr_start  = 0;
tr_end = 10*60*500;
    
subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
press_trunc = {'S1_';'S11_'}; % truncated pressure stimulus file

% define files to load
for sub = 1:length(subjects)
    cd(Pdata)
    subject = subjects(sub);
    subject = subject{:};

    s=[subject '*.vhdr'];
    eegfiles=dir(s);
    
% load data from each file  
for f=1:length(eegfiles)  

filename = eegfiles(f).name;
tempcnt = pop_loadbvCAB(Pdata, filename);

v2 = findstr(filename,'_2');
v3 = findstr(filename,'_3');

%correct bad channels

if strcmp(num2str(subject), 'S4_') == 1
    tempcnt.data(22,:) = mean(tempcnt.data([34 36],:),1);
    %tempcnt.data(23,:) = mean(tempcnt.data([35 37],:),1);
    tempcnt.data(3,:) = mean(tempcnt.data([33 47],:),1);
    %tempcnt.data(4,:) = mean(tempcnt.data([34 48],:),1);
    tempcnt.data(5,:) = mean(tempcnt.data([35 49 41 43],:),1);
    %tempcnt.data(19,:) = mean(tempcnt.data([63 64],:),1);
elseif strcmp(subject, 'S5_') == 1
    tempcnt.data(35,:) = mean(tempcnt.data([18 5],:),1);
    tempcnt.data(21,:) = mean(tempcnt.data([33 35],:),1);
    %tempcnt.data(4,:) = mean(tempcnt.data([34 48],:),1);
    tempcnt.data(22,:) = mean(tempcnt.data([34 36],:),1);
    tempcnt.data(37,:) = mean(tempcnt.data([19 7],:),1);
    tempcnt.data(62,:) = mean(tempcnt.data([17 61 39 40],:),1);
elseif strcmp(num2str(subject), 'S8_') == 1
    %tempcnt.data(6,:) = mean(tempcnt.data([36 50],:),1);
    tempcnt.data(22,:) = mean(tempcnt.data([34 36],:),1);
    %tempcnt.data(34,:) = mean(tempcnt.data([4 17],:),1);
elseif strcmp(subject, 'S9_') == 1
    tempcnt.data(17,:) = mean(tempcnt.data([33 34],:),1);
%    tempcnt.data(21,:) = mean(tempcnt.data([33 35],:),1);
%    tempcnt.data(1,:) = mean(tempcnt.data([61 53],:),1);
%    tempcnt.data(19,:) = mean(tempcnt.data([37 38],:),1);
%    tempcnt.data(22,:) = mean(tempcnt.data([34 36],:),1);
%    tempcnt.data(24,:) = mean(tempcnt.data([44 63],:),1);
%    tempcnt.data(43,:) = mean(tempcnt.data([23 27],:),1);
%    tempcnt.data(64,:) = mean(tempcnt.data([19 20 45 46],:),1);
%    tempcnt.data(37,:) = mean(tempcnt.data([23 45 7 64],:),1);	
elseif strcmp(subject, 'S10_') == 1
    tempcnt.data(2,:) = mean(tempcnt.data([61 54],:),1);
end

event = tempcnt.event;
stim = struct2cell(event);
chans = tempcnt.chanlocs;
chans= struct2cell(chans);
chans = squeeze(chans(1,:));
offset = squeeze(stim(1,:,:));
offset=[offset{:}];
offset_stim = offset(2:length(offset));
offset_isi = [];
for i = 1:length(offset_stim)-1
    offset_isi(i) = offset_stim(i+1) - offset_stim(i);
end
off = find(offset_isi>10000);
off_num = numel(off);
stimtype = squeeze(stim(5,:,:))';
ii=find(strcmp(stimtype,'S  1'));

total_data = [];
data = zeros(64,tr_end,2);

cd(Pdata)

if offset(ii(end))-tr_end+1 < 0
    tr_end = offset(ii(end));
end

data(:,:,1)=tempcnt.data(:, offset(ii(end))-tr_end+1:offset(ii(end)));
data(:,:,2)=tempcnt.data(:, offset(ii(end)):offset(ii(end))+tr_end-1);

for i=1:size(data,2)/(5*500)
    data1(:,:,i)=data(:, (5*500)*(i-1)+1:(5*500)*i,1);
    data2(:,:,i)=data(:, (5*500)*(i-1)+1:(5*500)*i,2);
end

total_data = cat(3, data1,data2);

clear data data1 data2

cd(Presults)

if v2; session = 2; elseif v3; session = 3; end
x=size(total_data);
eval(['save ' num2str(subject) 'data_matrix_dim_' num2str(session) ' x']);

if exist('total_data','var')
    Nelectrodes = x(1);
    Nsamples = x(2);
    Nevents = x(3);
    total_data = reshape(total_data,Nelectrodes,Nevents*Nsamples);
    total_data = total_data(ele,:,:);
    eval(['save ' num2str(subject) 'data_segments' num2str(session) '.mat'  ' total_data']);
    b= jader(total_data, 30);
    s=[num2str(subject) 'b' num2str(session) '.mat'];
    save(s, 'b');
end

clear tempcnt total_data total_data1 total_data2 v2 v3

end
end % end of subjects loop