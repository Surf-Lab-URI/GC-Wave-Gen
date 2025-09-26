clear all
close all

hRes = 2.6336d-4; %m/pixel
vRes = 2.7152d-4; %m/pixel
pivRes = 40d-6; %m/pixel
pivsurfRes = 102d-6; %m/pixel
lfvToPivScale = hRes/pivRes;

load('\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\Lfv\ExpLC1_dt25ms_1_Lfv_1002.mat')

imgLfv_LD = correctLfvLensDist(imgLfv);
imgLfv_fin = correctLfvCamAngle(imgLfv_LD);

lfvSurf = findSurfaceLfvFabMarc(imgLfv_fin, 30);

figure, imagesc(imgLfv_fin), colormap(bone), caxis([0 1000]), hold on, plot(lfvSurf.z_s,'r')

s = lfvSurf.z_s;

s1 = s*vRes; % surface in m
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(imgLfv_fin,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
s2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)
%%
% ssrvs_cut = ssrvs(end-lr-t2+2:end-lr+1);
% ssrvs_fin = (ssrvs_cut-ud);

% figure, imagesc(xi,yi,imgLfv_fin), colormap(bone), caxis([0 1000]), hold on, plot(xt,s2,'r')
% 
% figure, imagesc(PivFuse), colormap(bone), caxis([0 500])

%% common point between lfv and piv

X_lfv = 985; 
X_lfv1 = X_lfv * lfvToPivScale;
X_pivInLfv = round(X_lfv1-2255);
% X_piv = 2255;
% X_piv_m = X_piv * pivRes;

s2_pivSample = s2(X_pivInLfv:X_pivInLfv+3785-1);

exp_name = 'LC1_dt25ms_1';
path1 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv1\';
path2 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv2\';
path = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\pivsurf\';
image_pair_number = 1002;
num_of_digits = 4;
image_letter = 'a';
u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
fusedIm = fuseImagesLC(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter ,tform);

xpiv = 0:pivRes:3784*pivRes;
ypiv = 0:pivRes:2048*pivRes;
figure, imagesc(xpiv,ypiv,fusedIm.fused_im), colormap(bone), caxis([0 400]), hold on, plot(s2_pivSample-0.3,'r')

% fused_im = 
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
s3 = imSurfLC.z_s * pivsurfRes;


figure, plot(s2_pivSample-mean(s2_pivSample)), hold on, plot(s3(end:-1:1)-mean(s3),'r')

figure, imagesc(imSurfLC.pivSurf_resized), colormap(bone)


















