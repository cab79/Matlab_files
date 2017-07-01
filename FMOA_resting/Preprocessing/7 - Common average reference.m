clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA3','OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};


for n=1:length(subjects)

subject = subjects(n);
subject = char(subject);

fnames={
    '_total_data_ICA.mat';
  };


for x=1:length(fnames);

fname=char(fnames(x));
subject=char(subjects(n));
fname= [subject fname(1:length(fname)-4)];
load(fname);
mcr= mean(total_data_ICA([1:2 4:30 33:64],:));
for i=1:64
    total_data_ICA(i,:)=total_data_ICA(i,:)- mcr;
end
eval(['save ' fname '_ca.mat'  ' total_data_ICA'])

end

end


