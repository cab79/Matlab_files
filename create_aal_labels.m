function lab = create_aal_labels(S)
% S must contain:
    % S.aal_path
    % S.cimg : cluster img file name
    % S.folder : folder containing the img and for saving

V=spm_vol(fullfile(S.folder,[S.cimg '.nii']));
img = spm_read_vols(V);
    
[pth,nme,ext] = fileparts(S.aal_path);

% get labels and their indices
xml = xml2struct(fullfile(pth,[nme '.xml']));
aalstruct = vertcat(xml.atlas.data.label{:});
labels = vertcat(aalstruct(:).name);
index = vertcat(aalstruct(:).index);
index = cellfun(@str2double,{index(:).Text});

aal = load_nii(S.aal_path);
[ui,IA,IB] = unique(img(img>0));
lab = cell(length(ui),1);
for i = 1:length(ui)
    reg = reshape(img==ui(i),size(img));
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
        mni = cor2mni(round(cent{cc}([2 1 3])), V.mat); % convert to MNI
        lab{i,3}(cc,:) = mni; 
    end
    lab{i,4} = length(cent);
end
save(fullfile(S.folder,'aal_labels.mat'),'lab');