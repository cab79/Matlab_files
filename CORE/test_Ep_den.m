%test EP_Den


for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    if strcmp(C{2},'2')
        timebin = timebins(1,:);
    elseif strcmp(C{2},'4')
        timebin = timebins(2,:);
    end
    
    elec = [34 70];
    datin = reshape(squeeze(mean(EEG.data(elec,:,:),1)),EEG.trials*EEG.pnts,1);
    fid = fopen('testdat.asc','wt');
    fprintf(fid,'%2.12f\n',datin);
    fclose(fid);
    
    elec = [16];
    datin = reshape(squeeze(mean(EEG.data(elec,:,:),1)),EEG.trials*EEG.pnts,1);
    fid = fopen('testdat16.asc','wt');
    fprintf(fid,'%2.12f\n',datin);
    fclose(fid);
    
    
    handlesin.par.samples = EEG.pnts;
    handlesin.par.stim = -timebin(1)*Sr+1;
    handlesin.par.sr = Sr;
    handlesin.par.scales = 5;
    handlesin.par.max_trials = 100;
    handlesin.par.max_contour = 100;
    
    [datout H] = EP_den(datin,handlesin);
    
    clear dat
    
end
