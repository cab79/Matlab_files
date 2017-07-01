clear all
close all
addpath('F:\Dell\bloodA\Image_analysis_files\Matlab code\Create_parametric_images');
anapath = 'C:\Data\PET-LEP\PET-EEG_similarity_analysis';
savepath='C:\Data\PET-LEP\PET-EEG_similarity_analysis';
run('C:\Data\Matlab\Matlab_files\PET-LEP\EEG scripts\loadsubj.m');
grplist = [2,3,... % healthy
            5,6];   % patient

%% PET data info
Pbase = 'F:\Dell\bloodA\Image_analysis_files\examples';
PsubPET = 'PET';
scans = {1 2 'diff'};
no_scans = length(scans);

%% EEG data info
folder_name = 'GS_grp_20170222T120001';
Ebase = fullfile('C:\Data\PET-LEP\SPM_source_images',folder_name);
no_cond = 2;
epeaks_start = [1:14];

%% ROI data info
Rbase = 'SPM analysis\Brain_atlas_Hammers2003';
Rfile = dir(fullfile(Pbase,Rbase,['rHammers*SPM5.nii'])); % should be same size and voxel dimension as PET and EEG images.
Rdata = fullfile(Pbase, Rbase, Rfile.name);


mkdir(fullfile(anapath,folder_name));
cd(fullfile(anapath,folder_name));

subjects = subjlists(grplist);

%regions = {
    %'hipp',1:2;
    %'Lamyg',3;
    %'Ramyg',4;
    %'Amyg',3:4;
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
    %'ins', 20:21;
    %'acc',24:25;
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
    %'subcal', 78:79;
    %'presub',80:81;
    %'lingual',64:65;
    %'cuneus',66:67;
    %'lat_occ',22:23;
    %'strg',52:53;
    %'pcg',60:61;
    %'sup_par',62:63;
    %'inf_lat_par',32:33;
    %'caud',34:35;
    %'nuc_acc',36:37;
    %'puta',38:39;
    %'thal', 40:41;
    %'pall',42:43;
    %'cc',44;
    %'sub_nig',74:75;
    %'lat_vent',45:46;
    %'lat_vent_temp',47:48;
    %'third_vent',49;

    %};

Nsub = size(subjects, 1);

%load ROI file
[Rpth,Rnme,Rext] = fileparts(Rdata);
if strcmp(Rext,'.img')
    Rhdr = HDRread(fullfile(Rpth,Rnme),'ieee-le');
    RnXY = prod(Rhdr.dim(2:3));
    RnZ = Rhdr.dim(4);
    Rhdr.dim(5:8) = [1 1 1 1];
elseif strcmp(Rext,'.nii')
    Rnii = load_nii(Rdata);
end

reg_start = unique(Rnii.img(:));
reg_start(reg_start==0)=[];

epeaks = epeaks_start;
reg = reg_start;
Nreg = length(reg_start);
no_peaks = length(epeaks_start);
count=0;
while no_peaks>4 && Nreg>5 
    count=count+1
    
    % 3 PET cond (1,2,diff)
    % 2 EEG cond
    % 14 time points
    % 83 ROIs
    % multiple source models
    % after each run, find the most common EEG cond, PET cond, timebins and
    % ROIs showing the most effects, then narrow down the search.
    iii=0;
    results = [];
    for sc = 1:no_scans
        for c = 1:no_cond
             for p = 1:no_peaks
                 Nsub=0;
                 for s = 1:length(subjects)
                     for s2 = 1:size(subjects{s,1},1) 
                        Nsub=Nsub+1;

                        %list EEG filenames
                        tmp_nme = subjects{s,1}{s2,1};
                        file = dir(fullfile(Ebase,['rspm12_' tmp_nme '_' num2str(epeaks(p)) '_*' num2str(c) '.nii'])); % use resliced versions starting with 'r'
                        fnames{Nsub,1} = fullfile(Ebase,file.name);

                        %identify and list PET file names
                        Pfiles = dir(fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,['s8wrh*nnls_Vd_mag-310_rsl.img']));
                        Pfiles_sub = dir(fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,['s8wrh*nnls_Vd_mag-310_rsl_SUBTRACTION.img']));
                        Pfiles = [Pfiles;Pfiles_sub];

                        if ismember(grplist(s),[2 5])
                            scan_ana = subjects{s,1}{s2,3};
                        elseif ismember(grplist(s),[3 6])
                            scan1 = subjects{s,1}{s2,3};
                            scan2 = [scans{1:2}];
                            scan2(scan1) = [];
                            scan_ana = {scan1,scan2,scans{3}};
                            scan_order = [scan1,scan2];
                            Pfiles(1:2) = Pfiles(scan_order);
                        end

                        if length(Pfiles) ~= length(scan_ana); errordlg('number of images is incorrect'); end;

                        if sc<2
                            fnames{Nsub,2} = fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,Pfiles(sc).name);
                        else
                            if length(Pfiles)>1
                                fnames{Nsub,2} = fullfile(Pbase, subjects{s,1}{s2,2}, PsubPET,Pfiles(sc).name);
                            else
                                continue
                            end
                        end

                        % load PET file
                        [pth nme ext] = fileparts(fnames{Nsub,2});
                        sdate = nme(8:15);
                        Pnii = load_nii(fnames{Nsub,2});
                        if ~all(size(Pnii.img) == size(Rnii.img)); errordlg('PET-ROI images are not the same size'); end;

                        % load EEG file
                        [pth nme ext] = fileparts(fnames{Nsub,1});
                        Enii = load_nii(fnames{Nsub,1});
                        if ~all(size(Enii.img) == size(Rnii.img)); errordlg('EEG-ROI images are not the same size'); end;

                        for r = 1:Nreg
                            P=Pnii.img(ismember(Rnii.img,reg(r)));
                            E=Enii.img(ismember(Rnii.img,reg(r)));

                            %figure;scatter(tiedrank(P),tiedrank(E));lsline
                            %figure;scatter(P,E);lsline

                            results(Nsub,r,p,c,sc) = corr(P,E,'type','Spearman');
                            iii=iii+1
                            %results(s+2,(1+3*r-3+i)) = {Prs/R_size};
                        end
                    end
                end
             end
        end
    end

    findmax = 0.5;
    zs=[];
    rmax = [];
    pmax = [];
    for sc = 1:no_scans
        for c = 1:no_cond
             for p = 1:no_peaks
                 for r = 1:Nreg
                    zs(r,p,c,sc) = mean(results(:,r,p,c,sc))/std(results(:,r,p,c,sc));
                 end
             end
             M = zs(:,:,c,sc);
             [C,I] = sort(M(:),'descend');
             [I1,I2] = ind2sub(size(M),I(1:ceil((length(I)*findmax))));
             rmax(c,sc,:) = I1;
             pmax(c,sc,:) = I2;
        end
    end

    rFreq = [unique(rmax),histc(rmax(:),unique(rmax))];
    [rvalues, rorder] = sort(rFreq(:,2),'descend');
    rsorted = rFreq(rorder,:);
    rchosen = reg(rsorted(rsorted(:,2)>max(rsorted(:,2))/2,1));

    pFreq = [unique(pmax),histc(pmax(:),unique(pmax))];
    [pvalues, porder] = sort(pFreq(:,2),'descend');
    psorted = pFreq(porder,:);
    pchosen = epeaks(psorted(psorted(:,2)>max(psorted(:,2))/2,1));
    
    save(fullfile(savepath,['results_' datestr(datetime,30) '.mat']),'results','rchosen','pchosen');
    
    epeaks = sort(pchosen);
    reg = sort(rchosen);
    if Nreg==length(reg) && no_peaks==length(epeaks)
        break
    end
    Nreg = length(reg);
    no_peaks = length(epeaks);

end

close all
tsh=[];
tsp=[];
for sc = 1:no_scans
    for c = 1:no_cond
          for p = 1:no_peaks
            for r = 1:Nreg
                [tsh(r,p,c,sc) tsp(r,p,c,sc)] = ttest(results(:,r,p,c,sc),0,'Alpha',0.001);
            end
          end
        figure
        data=tsh(:,:,c,sc).*zs(:,:,c,sc);
        imagesc(data)
        xticks = linspace(1, size(data, 2), numel(pchosen));
        set(gca, 'XTick', xticks, 'XTickLabel', pchosen)
        yticks = linspace(1, size(data, 1), numel(rchosen));
        set(gca, 'YTick', yticks, 'YTickLabel', rchosen)
        xlabel('peaks'); ylabel('regions');
        title(['scan' num2str(sc) '-' 'cond' num2str(c)])
        temp=tsp(:,:,c,sc);
        temp=sort(temp(:))';
        temp=temp(1:5)
    end
end

