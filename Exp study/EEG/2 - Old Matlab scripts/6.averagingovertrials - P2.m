subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

for subject = 1:length(subjects)


d1 = dir([char(subjects(subject)) '_total_data_ica.mat']);
d2 = dir([char(subjects(subject)) ' _block_info.mat']);
load(d1.name); load(d2.name);

for i=1:size(total_data_ICA,3)
    total_data_ICA(:,:,i)=(detrend2(squeeze(total_data_ICA(:,:,i))',250,250))';
end

total_data_ICA = blcorrect4(total_data_ICA, 250);

block_size = max(events_mat(:,1));

for block = 1:block_size

    for trial = 1:6
    ii = find(events_mat(:,1)==block & events_mat(:,2)==trial & events_mat(:,3)==1);
avg = squeeze(mean(total_data_ICA(:,:,ii),3));

name=[char(subjects(subject)) '_' num2str(trial) '_avgP2'];
%name=[char(subjects(subject)) '_' num2str(block) '_' num2str(trial) '_avg']; %%% use if there is more than one block
save(name, 'avg'); clear avg;

    end
end
end