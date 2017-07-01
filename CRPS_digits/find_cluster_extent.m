% find extent of cluster in a specific dimention
dim = 3; % e.g. 3 is time if cluster is sensorspace
nii=load_nii('cluster124.nii');
ext = find(any(reshape(nii.img,[],size(nii.img,dim)),1));
first = ext(1)
last = ext(end)