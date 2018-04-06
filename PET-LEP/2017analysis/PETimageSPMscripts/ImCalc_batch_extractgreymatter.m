clear all

Pbase = 'C:\bloodA\Image_analysis_files\examples';
cd(Pbase)
Spaths = spm_select(Inf, 'dir', 'Select subject directories');
Nsub = size(Spaths, 1);

PsubPET = 'PET';
PsubMR = 'MR';

%for r = 1:Nreg
    for s = 1:Nsub

        Spath = deblank(Spaths(s, :));

        cd(fullfile(Spath, PsubMR))
        Mfile = dir('wc1*.img');
        Mdata = fullfile(Spath, PsubMR, Mfile.name);

        cd(fullfile(Spath, PsubPET))
        Pfiles = dir('wr*nnls_K1_mag-320_med_i3_d3_rsl.img');

        for i = 1:length(Pfiles)
            Pdata = fullfile(Spath, PsubPET, Pfiles(i).name);
            input = {Pdata;Mdata};

            [pth nme ext] = fileparts(Pfiles(i).name);

            Calc = spm_imcalc_ui(input,fullfile(Spath, PsubPET, [nme '_greymatter' ext]),'(i1.*(i2>0))');
        end

    end
%end
