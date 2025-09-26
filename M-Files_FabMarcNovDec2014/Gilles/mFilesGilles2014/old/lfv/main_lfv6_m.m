tic,
clear all
close all

pivResReal = 40d-6;
pivRes = 1;%40d-6; %m/pixel
iws = 16; % initial widow size

imnum = '1150';
lim = 4700;

load(['\\beo\data\ExpLC2_dt7ms_1\RawImages\Lfv\ExpLC2_dt7ms_1_Lfv_' imnum '.mat'])

%% correct lfv
imgLfv_t = correctLfv1(imgLfv);
imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);

%% find surface on lfv
lfvSurf = findSurfaceLfv(imgLfv_t_cut);
figure, imagesc(imgLfv_t_cut), colormap(bone), caxis([0 2000]), hold on, plot((lfvSurf.z_s),'r')

%% interpolote surface to resolution of piv
%%
vRes =  1.3543e-04/pivResReal; %pixel de piv/pixel de lfv;  %lfv resolution
hRes = 1.2912e-04/pivResReal; %pixel de piv/pixel de lfv
%%
s = lfvSurf.z_s;
s1 = s*vRes; % surface in m
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(imgLfv_t_cut,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
ss2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)

% figure, plot(ss2)
s2 = size(imgLfv_t,1)*hRes - ss2;
% figure, plot(s2)
s2 = s2 - mean(s2);

%% function compute orbitals
%%  calculate orbitals
%VELOCITY RECONSTRUCTION
% xx = 0:1:13100-1;
% s2 = sin(pi*xx/5000);
dx = pivRes;
x_s = 0 : pivRes * iws : (13100-1)*pivRes;
z_s = 0 : pivRes * iws : 2048; %vertical coordinate up to 20cm with 1mm resolution (to be adjusted to match PIV)
u = zeros(length(z_s),length(x_s)-1);
w = zeros(length(z_s)-1,length(x_s)-1);
f = fft(s2); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx; %wavenumber of mode
for j = 1:length(z_s)
    g = 0;
    h = 0;
    for i = 1:floor(length(s2)/2)%/10
        g = g-fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i))+mean(s2)/floor(length(s2))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
        h = h+fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i))-mean(s2)/floor(length(s2))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
    end
    u(j,:) = g;
    w(j,:) = h;
end
%% end compute orbitals
%% compute sPSa and sPSb
exp_name = 'LC2_dt7ms_1';
path1 = ['\\beo\data\Exp' exp_name '\RawImages\Piv1\'];
path2 = ['\\beo\data\Exp' exp_name '\RawImages\Piv2\'];
path = ['\\beo\data\Exp' exp_name '\RawImages\Pivsurf\'];
image_pair_number = str2double(imnum);
num_of_digits = 4;
image_letter = 'a';

u1 = [ 598 1577 990]'; v1 = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u1 v1],[x y]);
imSurfLCa = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% s3 = imSurfLC.z_s * pivsurfRes;
imgPivsurfa = fliplr(imSurfLCa.pivSurf_resized);
sPSa = imSurfLCa.z_s(end:-1:1).*pivRes;
%%
image_letter = 'b';
imSurfLCb = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% s3 = imSurfLC.z_s * pivsurfRes;
imgPivsurfb = fliplr(imSurfLCb.pivSurf_resized);
sPSb = imSurfLCb.z_s(end:-1:1).*pivRes;
% figure,imagesc(imgPivsurf), colormap(bone), hold on, plot(imSurfLC.z_s(end:-1:1), 'r')
% 
%% function interpolate orbitals from surface following grid to cartesian grid

imPSa = imSurfLCa.pivSurf_resized;
sPS = sPSa;
sPS = size(imPSa,1)*pivRes - sPSa;
%%
uu = u(2:end,floor(lim/iws):floor((lim+3784)/iws));
ww = w(2:end,floor(lim/iws):floor((lim+3784)/iws));
% s2s = s2(lim:lim+3784);
uui = nan(size(uu));
wwi = nan(size(ww));

% xii = ;
for col = 1:size(uu,2)-2
%     keyboard
    surfTemp = sPS((col-1)*iws+iws/2);
    xii = surfTemp:-pivRes*iws:0;
    xxt = 0:pivRes*iws:surfTemp;
    uui(1:length(xii),col) = spline(xii,uu(1:length(xii),col),xxt);
    wwi(1:length(xii),col) = spline(xii,ww(1:length(xii),col),xxt);
end
%% end function
figure, imagesc(uui),colorbar
figure, imagesc(wwi),colorbar

% figure, imagesc(u), colorbar, caxis([-0.2 0.2])
% figure, imagesc(w), colorbar, caxis([-0.2 0.2])
beg = 1;
st = 5;
figure, quiver(beg:st:size(uui,2), beg:st:size(uui,1), uui(beg:st:size(uui,1),beg:st:size(uui,2)), wwi(beg:st:size(uui,1),beg:st:size(uui,2)))
% 

toc
%%
u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
image_letter = 'a';
fusedImLCa = fuseImagesLC(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter,tform);
IM1 = fusedImLCa.fused_im;
image_letter = 'b';
fusedImLCb = fuseImagesLC(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter,tform);
IM2 = fusedImLCb.fused_im;
%%
mask1 = imSurfLCa.mask;
mask2 = imSurfLCb.mask;
delxOrb1 = uui(1:end-1,1:end-3);
delxOrb = flipud(-delxOrb1);
delzOrb1 = wwi(1:end-1,1:end-3);
delzOrb = flipud(delzOrb1);
IM1 = fliplr(IM1);
IM2 = fliplr(IM2);
mask1 = fliplr(mask1);
mask2 = fliplr(mask2);
compVel_orb = PIV_FAB6_LfvOrb1 (IM1, IM2, mask1, mask2, delxOrb, delzOrb);
compVel_orb_zero = PIV_FAB6_LfvOrb1 (IM1, IM2, mask1, mask2, zeros(size(delxOrb)), zeros(size(delxOrb)));
% numel(find(dcor<0.5))/numenanl(dcor)
% numel(find(isnan(dcor(:))))/numel(dcor(:))
dcor = compVel_orb.dcor;
MASK = compVel_orb.MASK;
dcor1 = 50 * ones(size(dcor)); % matrix of 50s
dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
figure, imagesc(dcor1), caxis([0 1]), colorbar

% figure, imagesc(dcor_zero)
% mask3 = MASK;
% mask3(isnan(MASK)) = 50;
% % figure, imagesc(mask3)
% mask3(mask3==1) = 0;
% figure, imagesc(mask3)
