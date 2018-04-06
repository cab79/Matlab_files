clear all

files = {'gamma1.mat';
'gamma2.mat';
'beta.mat'; 
'alpha.mat'; 
'theta.mat'; 
'delta.mat'}; 

for f = 1:length(files)
    
load(files{f});
[pth,fname,ext] = fileparts(files{f});
eval(['data = ' fname ';']) %dimensions are [ele,subject,scan,time]
FvPv = [];
for e = 1:size(data,1)
    edat = squeeze(data(e,:,:,:));
    sd = size(edat);
    % reshape such that time is the largest factor
    edat = reshape(edat,1,sd(1)*sd(2)*sd(3));
    % create design matrix 
    subjects = repmat([1:sd(1)],1,sd(2)*sd(3));
    scans = repmat(sort(repmat([1:sd(2)],1,sd(1)), 'ascend'),1,sd(3));
    time = sort(repmat([1:sd(3)],1,sd(1)*sd(2)), 'ascend');
    design = [edat' time' scans' subjects'];
    
    [P F] = RMAOV2_mod(design,0.05,0);
    
    FvPv(e,1,1:3) = F;
    FvPv(e,2,1:3) = P;
    
end

eval(['save ' fname '_ANOVAele_stats FvPv;'])

end
