n_per_grp = [13 13]; % participant numbers per group
n_cond_within = 5; % number of within-subjects conditions
n_cond_betw = 2; % number of between-subjects conditions

n_grps = length(n_per_grp);
tot_n = sum(n_per_grp);
tot_files = tot_n*n_cond_within*n_cond_betw;
factmat = ones(tot_files,4);


factmat(f,2) = 