files = dir('*');

for f = 1:length(files)
    try
        movefile(files(f).name,[files(f).name(1:4) '_' files(f).name(5:end)]);
    end
end