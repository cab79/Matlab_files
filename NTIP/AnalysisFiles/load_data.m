dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
S.loadpath = 'C:\Data\NTIP\Raw'; %where the raw input data is stored
S.savepath = 'C:\Data\NTIP\Raw'; %where the .set output files will go
S.fileext = '*.vhdr'; % file extension of input data
S.parts = [1];
S.combine=0; 

files = dir(fullfile(loadpath,S.fileext));
for f = 1:length(files)
    dataimport_bv(loadpath,files(f).name,parts,combine,savepath)
end