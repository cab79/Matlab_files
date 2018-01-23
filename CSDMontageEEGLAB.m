cd('C:\Data\Catastrophising study\Orig');
load chanlocs

trodes=[];
for site = 1:92
    trodes{site}=chanlocs(site).labels;
end;
trodes=trodes';

%% Get Montage for use with CSD Toolbox
Montage_92=ExtractMontage('10-5-System_Mastoids_EGI129.csd',trodes);
MapMontage(Montage_92);

%% Derive G and H!
[G,H] = GetGH(Montage_92);

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_92;                             % use generic variable name
save M:\Matlab\Matlab_files\CRPS_digits\CSDmontage_92.mat G H Montage; % save variables to Matlab file