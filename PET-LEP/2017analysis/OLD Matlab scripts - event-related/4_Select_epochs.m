clear all
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_';'P25_';'P35_'};
subject = 1; 
filename = 'total_data_ICA';
load ([subject filename '.mat']);
EEG = pop_loadset('S4_3.set','/scratch/cb802/Data/PET-LEP/Preprocessed data');
EEG.setname = subject;
EEG.filename = subject;
EEG.nbchan = 62;
eval(['EEG.data = ' filename ';']);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.setname = 'set';
EEG.nbchan = 62;
ALLEEG(1) = EEG;
CURRENTSET = 1;
pop_eegplot(ALLEEG,1,1,0);
reject = find(ALLEEG.reject.rejmanual == 1);