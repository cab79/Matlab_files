
addpath(genpath(pwd));

a=['root=''' pwd ''';'];
edit rootdir.m
FID=fopen('rootdir.m','a');
fprintf(FID, '%s', a);
fclose(FID);