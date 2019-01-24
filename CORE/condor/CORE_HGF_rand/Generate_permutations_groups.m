ngrps = [22,22];
nrand = 100;
save_dir = 'C:\Data\CORE\condor';

subs = [];
for g = 1:length(ngrps)
    subs = [subs g*ones(1,ngrps(g))];
end

urows=[];
i=0;
while size(urows,1)<nrand
    i=i+1;
    disp(num2str(i))
    grps_rand(i,:) = subs(randperm(length(subs)));
    urows = unique(grps_rand,'rows');
end

save(fullfile(save_dir,'grps_rand.mat'),'grps_rand');