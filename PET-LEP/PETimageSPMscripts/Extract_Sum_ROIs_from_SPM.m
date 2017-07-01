clear all
close all
Pbase = 'C:\bloodA\Image_analysis_files\examples';
PsubR = 'C:\bloodA\Image_analysis_files\examples\SPM analysis\SPM_ROIs';
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
    'h00625', 1;
    'h00627', 2;
    'h00629', 1;
    'h00642', 2;
    'h00644', 2;
    'h00645', 2;
    'h00646', 1;
    'h00650', 2;
    'h00658', 1;
    'h00660', 2;
    'h00663', 1;
    'h00667', 1;
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

regions = {
    %'caudate_McGsens'
    %'putamen_McGsens'
    %'S1_McGaff'
    %'mPFC_PCS'
    %'motor_PCS'
    %'ACC_maxpain'
    %'PFC_maxpain'
    %'Amyg_maxpain'
    %'motor_PCS_diff'
    %'PCC_P2'
    %'motor_pass'
    %'motor2_pass'
    %'PFC_pass'
    %'PCC_pass'
    %'temp_pass'
    %'caudate_painthresh'
    'ROI'
    };

filetypes = {'Vd', 'K1'};

Nsub = size(subjects, 1);
Nreg = size(regions, 1);


for f = 1%:length(filetypes)
    
    results = cell(Nsub+2,3*Nreg+1);
    results(3:end,1) = subjects(:,1);
    results(1,2:end) = reshape(repmat(regions(:,1)',3,1),1,3*Nreg);
    results(2,2:end) = reshape(repmat({'1','2','diff'},1,Nreg),1,3*Nreg);
    results(3:end,2:end) = num2cell(NaN(size(results(3:end,2:end))));
    
    for c = 1:size(results,2)
        results{1,c} = [filetypes{f} results{1,c} '_' results{2,c}];
    end

    for s = 1:Nsub

        Spath = fullfile(Pbase, subjects{s,1}, PsubPET);
        cd(Spath)
        Pfiles = dir(['s8wrh*nnls_' filetypes{f} '_mag-320*_rsl_MEANSUB.img']);

        for i = 1:length(Pfiles)

            scan1 = subjects{s,2};
            scan2 = scans;
            scan2(scan1) = [];
            scan_order = [scan1, scan2];

            Pdata = fullfile(Spath, Pfiles(scan_order(i)).name);
            [pth nme ext] = fileparts(Pdata);
            sdate = nme(8:15);
            hdr = HDRread(fullfile(pth,nme),'ieee-le');
            nXY = prod(hdr.dim(2:3));
            nZ = hdr.dim(4);
            hdr.dim(5:8) = [1 1 1 1];
            nframes = hdr.dim(5);


            for r = 1:size(regions,1)
                P = [];
                R = [];
                Prs=0;
                R_size = 0;

                cd(PsubR)
                Rfile = dir(['*' regions{r} '.img']);
                Rdata = fullfile(PsubR, Rfile.name);
                [Rpth,Rnme] = fileparts(Rdata);
                Rhdr = HDRread(fullfile(Rpth,Rnme),'ieee-le');
                RnXY = prod(Rhdr.dim(2:3));
                RnZ = Rhdr.dim(4);
                Rhdr.dim(5:8) = [1 1 1 1];

                if RnXY ~= nXY || RnZ ~= nZ; errordlg('images are not the same size'); end;


                for ipl = 1:nZ
                    P(:,:) = read_plane(Pdata,hdr,ipl,'ieee-le');
                    R(:,:) = read_plane(Rdata,Rhdr,ipl,'ieee-le');
                    P(isnan(P))=0;
                    Pr = sum(P .* (R>0));
                    R_size = R_size + length(find(R>0)==1);
                    Prs = Prs+Pr;
                end
                results(s+2,(1+3*r-3+i)) = {Prs/R_size};
                R_sizes(r) = R_size; 
            end

        end

    end

    for r = 1:Nreg
        for s = 1:Nsub
            results{s+2,(1+3*r-3+3)} = (results{s+2,(1+3*r-3+1)} - results{s+2,(1+3*r-3+2)})*100/results{s+2,(1+3*r-3+2)};
        end
    end

    eval(['results_' filetypes{f} ' = results;']);
    eval(['save ROIs_' filetypes{f} '.mat results_' filetypes{f} ';']);
end

