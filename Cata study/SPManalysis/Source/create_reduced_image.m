function reduced_img = create_reduced_image(S,clusname,folder,saveimg,ucol)

Cnii = load_nii(fullfile(folder,clusname));

img = Cnii.img;
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
    RCnii = Cnii;
    RCnii.img = reduced_img;
    save_nii(RCnii,fullfile(folder,'reduced.nii'))

    NCnii = Cnii;
    NCnii.img = new_clus_img;
    save_nii(NCnii,fullfile(folder,'new_clus.nii'))
    
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
        [ui,IA,IB] = unique(RCnii.img(RCnii.img>0));
        lab = cell(length(ui),1);
        for i = 1:length(ui)
            reg = reshape(RCnii.img==ui(i),size(RCnii.img));
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
        end
        save(fullfile(folder,'aal_labels.mat'),'lab');
    end

end