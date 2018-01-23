cd('C:\Data\Catastrophising study\Orig');
load chanlocs

trodes=[];
for site = 1:length(chanlocs)
    trodes{site}=chanlocs(site).labels;
end;
trodes=trodes';

%% Get Montage for use with CSD Toolbox
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129.csd',trodes);
MapMontage(Montage_64);

%% Derive G and H!
[G,H] = GetGH(Montage_64);

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save('C:\Data\Catastrophising study\Orig\CSDmontage_64.mat', 'G','H','Montage'); % save variables to Matlab file