%clear all
%close all
Pbase = 'C:\bloodA\Image_analysis_files\examples';
PsubR = 'SPM analysis\Brain_atlas_Hammers2003';
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
    %'hipp',1:2;
    %'Lamyg',3;
    %'Ramyg',4;
    'Amyg',3:4;
    %'ant_temp_med',5:6;
    %'ant_temp_lat',7:8;
    %'parahipp',9:10;
    %'sup_temp_post',11:12;
    %'mid_inf_temp',13:14;
    %'fusi',15:16;
    %'post_temp',30:31;
    %'sup_temp_ant',82:83;
    %'cereb',17:18;
    %'brainstem' 19;
    %'Lins', 20;
    %'Rins', 21;
    'ins', 20:21;
    'acc',24:25;
    %'pcc',26:27;
    %'mid_front',28:29;
    %'precentral',51:52;
    %'ant_orb',54:55;
    %'inf_front',56:57;
    %'sup_front',58:59;
    %'med_orb',68:69;
    %'lat_orb',70:71;
    %'post_orb',72:73;
    %'subg', 76:77;
    'subcal', 78:79;
    %'presub',80:81;
    %'lingual',64:65;
    %'cuneus',66:67;
    %'lat_occ',22:23;
    %'strg',52:53;
    %'pcg',60:61;
    'sup_par',62:63;
    'inf_lat_par',32:33;
    'caud',34:35;
    'nuc_acc',36:37;
    'puta',38:39;
    'thal', 40:41;
    %'pall',42:43;
    %'cc',44;
    %'sub_nig',74:75;
    %'lat_vent',45:46;
    %'lat_vent_temp',47:48;
    %'third_vent',49;
    
    
    
    
    };

Nsub = size(subjects, 1);
Nreg = size(regions, 1);

results = cell(Nsub+2,3*Nreg+1);
results(3:end,1) = subjects(:,1);
results(1,2:end) = reshape(repmat(regions(:,1)',3,1),1,3*Nreg);
results(2,2:end) = reshape(repmat({'1','2','diff'},1,Nreg),1,3*Nreg);
results(3:end,2:end) = num2cell(NaN(size(results(3:end,2:end))));

for c = 1:size(results,2)
    results{1,c} = ['Vd_' results{1,c} '_' results{2,c}];
end

cd(fullfile(Pbase, PsubR))
Rfile = dir(['Hammers*SPM5_rsl.img']);

Rdata = fullfile(Pbase, PsubR, Rfile.name);
[Rpth,Rnme] = fileparts(Rdata);
Rhdr = HDRread(fullfile(Rpth,Rnme),'ieee-le');
RnXY = prod(Rhdr.dim(2:3));
RnZ = Rhdr.dim(4);
Rhdr.dim(5:8) = [1 1 1 1];

for s = 1:Nsub
	
    Spath = fullfile(Pbase, subjects{s,1}, PsubPET);
    cd(Spath)
    Pfiles = dir(['s8wrh*nnls_Vd_mag-310_rsl.img']);
    
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
        
        if RnXY ~= nXY || RnZ ~= nZ; errordlg('images are not the same size'); end;
          
        for r = 1:size(regions,1)
            P = [];
            R = [];
            Prs=0;
            R_size = 0;
            for ipl = 1:nZ
                P(:,:) = read_plane(Pdata,hdr,ipl,'ieee-le');
                R(:,:) = read_plane(Rdata,Rhdr,ipl,'ieee-le');
                P(isnan(P))=0;
                Pr = sum(P .* ismember(R,regions{r,2}));
                R_size = R_size + length(find(ismember(R,regions{r,2})==1));
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

results{Nsub+3,1} = 'corrcoef';
for r = 1:Nreg*3
    x = cell2mat(results(3:Nsub+2,r+1));
    y = pred;
    z = covs;
    nans = find(isnan(y)==1);
    nonans = 1:length(y);
    nonans(nans) = [];
    x = x(1:length(y));
    x = x(nonans);
    y = y(nonans);
    z = z(nonans,:);
    
    if ~isempty(covs)
        [b,bint,r1,rint,stats] = regress(x,z(:,1));
        [b,bint,r2,rint,stats] = regress(r1,z(:,2));
        [b,bint,r3,rint,stats] = regress(r2,z(:,3));
        [b,bint,r4,rint,stats] = regress(r3,z(:,4));
        [b,bint,r5,rint,stats] = regress(r4,z(:,5));
        %[b,bint,r6,rint,stats] = regress(r5,Resid(nonans,6));
    else 
        r5 = x;
    end
    
    corr = corrcoef(y,r5);
    results{Nsub+3,r+1} = corr(1,2);
end

