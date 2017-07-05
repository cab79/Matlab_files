% directory containing all image directories
img_dir = 'C:\Data\Catastrophising study\SPMdata\sensorimages';
% normlisation expression
expr = 'log(i1+sqrt(i1.^2+1))'; % arsinh (https://en.wikipedia.org/wiki/Inverse_hyperbolic_function) to reduce outliers
% delete original image?
delete_orig = 0;
% overwrite if there is an image with the same filename/path?
overwrite = 0;

subdir = dir(fullfile(img_dir,'t-*'));
for d = 1:length(subdir)
    imgs = dir(fullfile(img_dir,subdir(d).name,'scond*'));
    for i = 1:length(imgs)
        img = fullfile(img_dir,subdir(d).name,imgs(i).name);
        imgo = fullfile(img_dir,subdir(d).name,['n' imgs(i).name]);
        if overwrite || ~exist(imgo,'file')
            Output = spm_imcalc_ui(img,imgo,expr);
        end
        if delete_orig==1
            delete(img);
        end
    end
end
        
        