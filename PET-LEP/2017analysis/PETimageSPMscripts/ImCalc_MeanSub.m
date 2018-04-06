Pbase = 'C:\bloodA\Image_analysis_files\examples';
PsubPET = 'PET';
scans = [1 2];

% 1 = pain scan first, 2 = pain scan second
subjects = {
    'h00536', 1;
    'h00569', 2;
    'h00593', 1;
    'h00596', 1;
    'h00601', 1;
    'h00606', 1;
    'h00617', 2;
    'h00627', 2;
    'h00629', 1;
    'h00642', 2;
    'h00644', 2;
    'h00645', 2;
    'h00646', 1;
    'h00650', 2;
    'h00660', 2;
    'h00663', 1;
    'h00674', 2;
    'h00344', 2;
    'h00404', 1;
    'h00405', 2;
    'h00406', 1;
    'h00432', 2;
    'h00437', 2;
    'h00440', 1;
    'h00566', 1;
    'h00576', 2;
    };

Nsub = size(subjects, 1);

input ={};

for s = 1:Nsub    
    Spath = fullfile(Pbase, subjects{s,1}, PsubPET);
    cd(Spath);
    Pfiles = dir(['s12wrh*nnls_Vd_mag-320_med_i3_d3_rsl.img']);
    Mfile = dir(['*_MEAN.img']);
    if isempty(Pfiles); errordlg([subjects{s,1} ' empty']);end
    scan1 = subjects{s,2};
    scan2 = scans;
    scan2(scan1) = [];
    scan_order = [scan1, scan2];
    for i = 1:size(Pfiles,1)
        Pdata = fullfile(Spath, Pfiles(scan_order(i)).name);
        input{1,1} = Pdata;
    
        Mdata = fullfile(Spath, Mfile.name);
        input{2,1} = Mdata;

        [pth nme ext] = fileparts(Pdata);
        Calc = spm_imcalc_ui(input,fullfile(Spath, [nme '_MEANSUB.img']),'i1-i2');
        
    end
    clear Pfiles input;
end


