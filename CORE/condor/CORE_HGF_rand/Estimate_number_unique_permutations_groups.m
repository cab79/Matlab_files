grps = [ones(1,22) 2*ones(1,22)];

for i = 1:100000
    disp(num2str(i))
    grps_rand(i,:) = grps(randperm(length(grps)));
end

urows = unique(grps_rand,'rows');

size(urows,1)