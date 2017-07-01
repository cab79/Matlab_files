function Create_pearsonsR_matrices(filepath,basename,conntype, trunc)
fprintf('Processing correlation matrix %s',basename);
if isunix
    filepath = '/scratch/cb802/Data/CIP_resting/fMRI';
else
    filepath = 'W:\Data\CIP_resting\fMRI';
end
nsurr=50;
[num,txt,raw] = xlsread('ROIs.xlsx');
ROInames = txt(3:162,7);

matrix2 = zeros(length(ROInames),length(ROInames));

for r = 1:length(ROInames)
    data(:,r) = load([basename ROInames{r} '.txt']);
end
%surrdata = phaseran(data,nsurr);

if trunc>0
    data = data(1:trunc,:);
end

matrix = corr(data);

matrix2(:) = 0;
matrix2 = tril(matrix,-1)+tril(matrix,-1)';
matrix = matrix2;

bootmat=[];

save(fullfile(filepath, conntype, [basename conntype '.mat']),'matrix','bootmat');

