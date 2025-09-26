% Function findSurfaceLC
%% Object:
% resizes and correct the angle of pivSurf LIF image to match fused PIV image, finds surface on resized pivSurf 
% and generates a binary mask for PIV algo
%% Arguments: 
% experiment name, path of data, image pair number, number of digits in
% image pair number, image letter, fused image, tformSurf
%% Result: 
% experiment number, image pair number, pivSurf resized, x-coords of surface, y-coords of
% surface, binary mask (1s and Nans)
%% Author:
% Marc Buckley, Gilles Bouille
%% Last update:
% 08/16/2013
%% Example:
% exp_name = '3';
% path = ['\\Afsx1\piv2\Exp' exp_name '\' 'test_ims\'];
% image_pair_number = 440;
% num_of_digits = 4;
% image_letter = 'a';
% u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
% fusedIm = fuseImages(exp_name, path, image_pair_number, num_of_digits, image_letter);
% imSurf = findSurface(exp_name, path, image_pair_number, num_of_digits, image_letter, fusedIm.fused_im,tformSurf);
% ans = 
%            exp_name: '5'
%         im_pair_num: 2130
%           im_letter: 'a'
%     pivSurf_resized: [2040x3944 double]
%                 x_s: [1x3944 double]
%                 z_s: [1x3944 double]
%                mask: [2040x3944 double]
%%
function imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf)
%[imSurf.exp_name, imSurf.im_pair_num, imSurf.im_letter, imSurf.pivSurf_resized, imSurf.x_s, imSurf.z_s, imSurf.mask]

pivSurf_struc = load([path 'Exp' exp_name '_PivSurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
pivSurf = fliplr(pivSurf_struc.imgPivsurf);
%% Correction (angle et taille) de Pivsurf
imgPivsurf_t=imtransform(pivSurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);

%% Detection de la surface
img = imgPivsurf_t;
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
%locate outliers and nan them
surfstd = std(surface);
surfmean = mean(surface);
surface(abs(surface-surfmean)>3*surfstd)=nan;
smth_surface = smoothn(surface,1000);
%
for i=1:size(img,2)
    mask(round(smth_surface(i)):2048,i) = 1;
end
% 
imSurfLC.exp_name = exp_name;
imSurfLC.im_pair_num = image_pair_number;
imSurfLC.im_letter = image_letter;
imSurfLC.pivSurf_resized = imgPivsurf_t;
imSurfLC.x_s = 1:3785;
imSurfLC.z_s = smth_surface;
imSurfLC.mask = mask;