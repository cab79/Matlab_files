d = 'Q:\CORE\EEG\ana\prep\cleaned_tm';

%newdir = {'epoched', 'cleaned_tm', 'merged_cleaned', 'manrej', '1st_ICA', '2nd_ICA', 'ALLEEG', 'merged', 'cleaned'}
newdir = {'_2_', '_4_'};

for n = 1:length(newdir)
    files = dir(fullfile(d,['*' newdir{n} '*']));
    for f = 1:length(files)
        if ~exist(fullfile(d,newdir{n}),'dir')
            mkdir(fullfile(d,newdir{n}))
        end
        movefile(fullfile(d,files(f).name),fullfile(d,newdir{n},files(f).name));
    end
end