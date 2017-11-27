function reduced_img = create_reduced_image(S,clusname,folder,saveimg,ucol)

Hnii=spm_vol(fullfile(folder,clusname));
img = spm_read_vols(Hnii);

if isempty(ucol)
    [~,clusfield,~] = fileparts(clusname);
    ucol = S.wf.(clusfield).ucol;
end

[uimg,iU,iI] = unique(img);
reduced_img = zeros(size(img));
new_clus_img = zeros(size(img));
for u = 1:length(uimg)
    if uimg(u)==0
        continue
    else
        reduced_img(iI==u) = ucol(uimg(u));
        % when does ucol have multiple repetitions of it's value?
        if sum(ucol==ucol(uimg(u)))>1
            new_clus_img(iI==u) = ucol(uimg(u));
        end
    end
end
if saveimg
    RCnii=Hnii;
    RCnii.fname=fullfile(folder,'reduced.nii');
    spm_write_vol(RCnii,reduced_img);

    NCnii=Hnii;
    NCnii.fname=fullfile(folder,'new_clus.nii');
    spm_write_vol(NCnii,new_clus_img);
    
    % name by AAL regions
    if S.use_aal
        [pth,nme,ext] = fileparts(S.aal_path);

        % get labels and their indices
        xml = xml2struct(fullfile(pth,[nme '.xml']));
        aalstruct = vertcat(xml.atlas.data.label{:});
        labels = vertcat(aalstruct(:).name);
        index = vertcat(aalstruct(:).index);
        index = cellfun(@str2double,{index(:).Text});

        aal = load_nii(S.aal_path);
        [ui,IA,IB] = unique(reduced_img(reduced_img>0));
        lab = cell(length(ui),1);
        for i = 1:length(ui)
            reg = reshape(reduced_img==ui(i),size(reduced_img));
            regaal = reg.*double(aal.img);
            lab{i,1} = unique(regaal(regaal>0));
            temp = {};
            for j = 1:length(lab{i,1})
                temp{1,j} = labels(find(index==lab{i,1}(j))).Text;
            end
            lab{i,2}=strjoin(temp,', ');
            if isempty(lab{i,2})
                lab{i,2} = 'unknown';
            end
            % get centroid coordinates
            STATS = regionprops(reg,'Centroid');
            cent = {STATS.Centroid};
            for cc = 1:length(cent)
                mni = cor2mni(round(cent{cc}([2 1 3])), RCnii.mat); % convert to MNI
                lab{i,3}(cc,:) = mni; 
            end
            lab{i,4} = length(cent);
        end
        save(fullfile(folder,'aal_labels.mat'),'lab');
    end

end