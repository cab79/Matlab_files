function test_audiospeed

cd('C:\Matlab_files\NTIP\SCIn');
outfile = fopen('latencytest.txt', 'w');
soundfilename = 'audioplayer';
load temp

for i = 1:10
    i
    t1=GetSecs;
    play(temp);
	t2=GetSecs;
	latency=t2-t1;
	fprintf(outfile, '%s,%2.4f\n', soundfilename, latency); 
    pause(1)
end

fclose(outfile);