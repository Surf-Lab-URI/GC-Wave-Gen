% Function findSurface
%% Object:
% resizes pivSurf LIF image to match fused PIV image, finds surface on resized pivSurf 
% and generates a binary mask for PIV algo
%% Arguments: 
% experiment name, path of data, image pair number, number of digits in
% image pair number, image letter, fused image, piv resolution, pivSurf
% resolution
%% Result: 
% experiment number, image pair number, pivSurf resized, x-coords of surface, y-coords of
% surface, binary mask (1s and Nans)
%% Author:
% Marc Buckley
%% Last update:
% 08/08/2013
%% Example:
% exp_name = '3';
% path = ['\\Afsx1\piv2\Exp' exp_name '\' 'test_ims\'];
% image_pair_number = 440;
% num_of_digits = 4;
% image_letter = 'a';
% fusedIm = fuseImages(exp_name, path, image_pair_number, num_of_digits, image_letter);
% imSurf = findSurface(exp_name, path, image_pair_number, num_of_digits, image_letter, fusedIm.fused_im);
% ans = 
%            exp_name: '5'
%         im_pair_num: 2130
%           im_letter: 'a'
%     pivSurf_resized: [2040x3944 double]
%                 x_s: [1x3944 double]
%                 z_s: [1x3944 double]
%                mask: [2040x3944 double]
%%
function imSurf = findSurface(exp_name, path, image_pair_number, num_of_digits, image_letter, fused_im, piv_res, pivsurf_res)
%[imSurf.exp_name, imSurf.im_pair_num, imSurf.im_letter, imSurf.pivSurf_resized, imSurf.x_s, imSurf.z_s, imSurf.mask]
%% resize pivsurf image to match fuse image
bool = false;
surfScalingFac = pivsurf_res/piv_res;
pivSurf_struc = load([path 'Exp' exp_name '_PivSurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
pivSurf = pivSurf_struc.imgPivsurf;
pivSurf1 = imresize(pivSurf,surfScalingFac,'Method', 'nearest');
[t1,t2] = size(fused_im);

% pivsurf is bigger than pivfused
lr = 255;
ud = 1052-5;

% take fraction of pivsurf that matches fused_piv
% pivSurf2 = fliplr(pivSurf1);
% img1 = pivSurf1(ud:ud+t1-1,lr:lr+t2-1);
pivSurf2 = pivSurf1(ud:ud+t1-1,end-lr-t2+2:end-lr+1);
img = pivSurf2;
%
smth_vert = nan(size(img));
grad_vert = nan(size(img));
mask = nan(size(img));
surface = nan(1,size(img,2));
% 
for i=1:size(img,2)
    imgi = img(:,i);
    %locate outliers and nan them
    imgistd = std(imgi);
    imgimean = mean(imgi);
    imgi(abs(imgi-imgimean)>3*imgistd)=nan;
    %smooth each column
    smth_vert(:,i) = smoothn(imgi,10000);
    %compute gradient on each column
    grad_vert(:,i) = gradient(smth_vert(:,i));
    gv = grad_vert(:,i);
    [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/2, 'npeaks',1);
    surface(i) = locs;
    clear imgi gv locs pks
end
%
%locate outliers and nan them
surfstd = std(surface);
surfmean = mean(surface);
surface(abs(surface-surfmean)>3*surfstd)=nan;
smth_surface = smoothn(surface);
%
for i=1:size(img,2)
    mask(1:round(smth_surface(i)),i) = 1;
end
% 
imSurf.exp_name = exp_name;
imSurf.im_pair_num = image_pair_number;
imSurf.im_letter = image_letter;
imSurf.pivSurf_resized = pivSurf2;
imSurf.x_s = 1:t2;
imSurf.z_s = smth_surface;
imSurf.mask = mask;