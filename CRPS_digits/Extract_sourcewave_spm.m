clear all
close all
dname = pwd;

cd(dname);
grplist = [39 40 41 42]; sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
%grplist = [33 34 31 32]; sublist_side = {'L','R','L','R'}; %Affected vs %unaffected exp2
cond_nme = {'1','2','3','4','5'};
hand_nme = {'L','R'};
no_cond = length(cond_nme); % no of conditions per data file (arm)

%Rdata = '/scratch/cb802/Brain_atlas_Hammers2003/Hammers_mith_atlas_n30r83_SPM5_rsl_'';
Rdata = 'W:\Brain_atlas_Hammers2003\';
%Rdata = 'W:\Brain_atlas_Hammers2003\';



%Sd = [Sdata '.hdr'];
    
%regions = {%'L_postcen',60;
%	'R_postcen',61;%
	%'L_sup_par',62;
	%'R_sup_par',63;
	%'L_inf_par',32;
%	'R_inf_par',33;
    %'peak4clus',[];
%    };

regions = {
    'Rpar',[];
    'Lpar',[];
    'Rfron',[];
    'Lfron',[];
    'R_postcen',[];
   'L_postcen',[];
    'S1',[];
    };

%regions = {
%    'MSP_Rinfpar',[];
%    'MSP_Rsuppar',[];
%    'MSP_Lsuppar',[];
%    };

Nreg = size(regions,1);

loadsubj
subjects = subjlists(grplist);

Ns=0;
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        Ns=Ns+1;
        tmp_nme = subjects{s,1}{s2,1};
        tmp_nme = strrep(tmp_nme, '.left', '_left');
        tmp_nme = strrep(tmp_nme, '.Left', '_left');
        tmp_nme = strrep(tmp_nme, '.right', '_right');
        tmp_nme = strrep(tmp_nme, '.Right', '_right');
        tmp_nme = strrep(tmp_nme, '.flip', '_flip');
        tmp_nme = strrep(tmp_nme, '.Flip', '_flip');
        tmp_nme = strrep(tmp_nme, '.aff', '_aff');
        tmp_nme = strrep(tmp_nme, '.Aff', '_aff');
        tmp_nme = strrep(tmp_nme, '.Unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '.unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '_Left', '_left');
        tmp_nme = strrep(tmp_nme, '_Right', '_right');
        tmp_nme = strrep(tmp_nme, '_Flip', '_flip');
        tmp_nme = strrep(tmp_nme, '_Aff', '_aff');
        tmp_nme = strrep(tmp_nme, '_Unaff', '_unaff');
        tmp_nme = strrep(tmp_nme, '.Exp1', '_Exp1');
        fnames{Ns,1} = ['spm8_' tmp_nme];
    end
end

for r = 1:Nreg
    %roi_o = maroi([Rdata regions{r,1} '.hdr']);
    %sp = native_space(roi_o);  %gives you the space of the image where the roi was defined
    %Pos = voxpts(roi_o,sp);
    %Posmm = sp.mat(1:3,1:3)*Pos + repmat(sp.mat(1:3,4),1,size(Pos,2));
    
    Rd = [Rdata regions{r,1} '.hdr'];
    
    V = spm_vol(Rd);
    maskdata = spm_read_vols(V);
    [x,y,z] = ind2sub(size(maskdata),find(maskdata));
    xyz = [x y z]';
    mni = vox2mni(V.mat,xyz);
    
    %Vs = spm_vol(Sd);
    %maskdataS = spm_read_vols(Vs);
    %[xS,yS,zS] = ind2sub(size(maskdataS),find(maskdataS));
    %xyzS = [xS yS zS]';
    %mniS = vox2mni(Vs.mat,xyzS);
    
    %mniCom = mni(:,ismember(mni',mniS','rows'));
    
    for i = 1:Ns
        fname = fnames{i};
        [pth nme ext] = fileparts(fname);
        D     = spm_eeg_load(fullfile(pwd,fname));
        
        %D = struct(D);
        %size(D.trials,1);
        %D.trials.labels
        
        D.val = 3;
        
        % Requires:
        %
         D.inv{D.val}.source.XYZ  = mni';  %- (n x 3) matrix of MNI coordinates
        %
        % Optional:
        %
         D.inv{D.val}.source.rad = 0; % - radius (mm) of VOIs (default 5 mm)
         D.inv{D.val}.source.label = repmat({regions{r,1}},size(mni,2),1); %- label(s) for sources (cell array)
         D.inv{D.val}.source.fname = [nme '_' regions{r,1}]; %- output file name
         D.inv{D.val}.source.type = 'evoked'; % - output type ('evoked'/'trials')

        [Ds, D] = spm_eeg_inv_extract_CAB(D);
        
        


    end
    
end





        