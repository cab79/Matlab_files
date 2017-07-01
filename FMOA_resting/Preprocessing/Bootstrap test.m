% INPUT Args:
subjects = {'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26'}%;'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33'};
NSub = length(subjects);
ele = [1:2 4:30 33:64];
conditions = 2;
data1_1 = [(1:NSub)' zeros(NSub,length(ele))];
data1_2 = [(1:NSub)' zeros(NSub,length(ele))];
data2_1 = [(1:NSub)' zeros(NSub,length(ele))];
data2_2 = [(1:NSub)' zeros(NSub,length(ele))];

%%%%%%%%%%
fnames={'_1_l_av.dat'
    '_2_l_av.dat'
    '_3_l_av.dat'
    '_4_l_av.dat'
  };

for n = 1:NSub
sub = subjects(n);
c1 = load([char(sub) char(fnames(1))]);
c2 = load([char(sub) char(fnames(2))]);
c3 = load([char(sub) char(fnames(3))]);
c4 = load([char(sub) char(fnames(4))]);

data1_1(n,2:length(ele)+1) = c1(:,[1:2 4:30 33:64]);
data1_2(n,2:length(ele)+1) = c2(:,[1:2 4:30 33:64]);
data2_1(n,2:length(ele)+1) = c3(:,[1:2 4:30 33:64]);
data2_2(n,2:length(ele)+1) = c4(:,[1:2 4:30 33:64]);

end

iterations = 1000; % the number of iterations to run
shuffles = [];   % result of call to randperm2. (if [], computes it as above)
doWilcox = 1;         % pick 1 to do non-parametric test instead of ttest
doSig = 1;          % whether to return the significance of theactual data as well

[p_boot,p_sig,shuffles] = bootstrap(iterations, data1_2, data2_2, shuffles, doWilcox, doSig);