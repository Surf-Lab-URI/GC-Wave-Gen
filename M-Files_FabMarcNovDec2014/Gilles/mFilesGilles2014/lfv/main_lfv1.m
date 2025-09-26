clear all
close all

% hRes = 2.6336d-4; %m/pixel
% vRes = 2.7152d-4; %m/pixel
pivRes = 40d-6; %m/pixel
pivsurfRes = 40d-6; %m/pixel
% lfvToPivScale = hRes/pivRes;

load('\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\Lfv\ExpLC1_dt25ms_1_Lfv_1002.mat')

imgLfv_t = correctLfv1(imgLfv);
imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
lfvSurf = findSurfaceLfvFabMarc(imgLfv_t_cut, 50);
figure, imagesc(imgLfv_t_cut), colormap(bone), caxis([0 2000]), hold on, plot(smoothn(lfvSurf.z_s),'r')

%%
vRes =  1.3543e-04; %m/pix;
hRes = 1.2912e-04; %m/pix;
%%

s = lfvSurf.z_s;

s1 = s*vRes; % surface in m
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(imgLfv_t_cut,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
s2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)
%%


%% common point between lfv and piv

% X_lfv = 985; 
% X_lfv1 = X_lfv * lfvToPivScale;
% X_pivInLfv = round(X_lfv1-2255);
% % X_piv = 2255;
% % X_piv_m = X_piv * pivRes;
% 
% s2_pivSample = s2(X_pivInLfv:X_pivInLfv+3785-1);

exp_name = 'LC1_dt25ms_1';
path1 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv1\';
path2 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv2\';
path = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\pivsurf\';
image_pair_number = 1002;
num_of_digits = 4;
image_letter = 'a';
% u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
% fusedIm = fuseImagesLC(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter ,tform);

% xpiv = 0:pivRes:3784*pivRes;
% ypiv = 0:pivRes:2048*pivRes;
% figure, imagesc(xpiv,ypiv,fusedIm.fused_im), colormap(bone), caxis([0 400]), hold on, plot(s2_pivSample-0.3,'r')

% fused_im = 
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% s3 = imSurfLC.z_s * pivsurfRes;
imgPivsurf =fliplr(imSurfLC.pivSurf_resized);
sPS = imSurfLC.z_s(end:-1:1).*pivsurfRes;
figure,imagesc(imgPivsurf), colormap(bone), hold on, plot(sPS, 'r')

sLfv = s2;
% figure, plot(s2_pivSample-mean(s2_pivSample)), hold on, plot(s3(end:-1:1)-mean(s3),'r')
figure, plot(sLfv-mean(sLfv))
figure, plot(sPS-mean(sPS),'r')

sPSStart = 4900;
figure, plot(sLfv-mean(sLfv)), hold on, plot(sPSStart:sPSStart+3784,sPS-mean(sPS),'r')

% figure, imagesc(imSurfLC.pivSurf_resized), colormap(bone)


















