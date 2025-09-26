tic,
clear all
close all

pivRes = 40d-6; %m/pixel
iws = 32; % initial widow size

imnum = '1500';

load(['\\beo\data\ExpLC2_dt7ms_1\RawImages\Lfv\ExpLC2_dt7ms_1_Lfv_' imnum '.mat'])

%% correct lfv
imgLfv_t = correctLfv1(imgLfv);
imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);

%% find surface on lfv
lfvSurf = findSurfaceLfv(imgLfv_t_cut);
figure, imagesc(imgLfv_t_cut), colormap(bone), caxis([0 2000]), hold on, plot((lfvSurf.z_s),'r')

%% interpolote surface to rsolution of piv
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

%%  calculate orbitals
%VELOCITY RECONSTRUCTION
% xx = 0:1:13100-1;
% s2 = sin(pi*xx/5000);
dx = pivRes;
x_s = 0 : pivRes * iws : (13100-1)*pivRes;
z_s = 0 : pivRes * iws : 20d-2; %vertical coordinate up to 20cm with 1mm resolution (to be adjusted to match PIV)
u=zeros(length(z_s),length(x_s)-1);w=zeros(length(z_s)-1,length(x_s)-1);
f=fft(s2); %Fourrier Modal decomposition
fa=2*abs(f)/length(f); %amplitude of mode
fp=angle(f); %phase of mode
k=[0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx; %wavenumber of mode
for j=1:length(z_s)
g=0;
h=0;
for i=1:floor(length(s2)/2)%/10
g=g-fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i))+mean(s2)/floor(length(s2))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
h=h+fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i))-mean(s2)/floor(length(s2))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
end
u(j,:)=g;
w(j,:)=h;
end

figure,imagesc(u), colorbar, caxis([-0.2 0.2])
figure,imagesc(w), colorbar,caxis([-0.2 0.2])
beg = 2;
figure, quiver(beg:10:409,beg:5:157,u(beg:5:157,beg:10:409),w(beg:5:157,beg:10:409))

%%
%Surface RECONSTRUCTION at t
x_s2=x_s;%0:13100;
dx2=dx;
t=7d-3;
f=fft(s2);
fa=2*abs(f)/length(f);
fp=angle(f);
k=[0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx2;
omega=sqrt(9.81*k);
g=0;
for i=1:floor(length(s2)/2)
g=g+(fa(i))*cos(k(i)*x_s2(1:end-1)+fp(i)-omega(i)*t)-mean(s2)/floor(length(s2)/2);
end
g=interp1(x_s2(1:end-1),g,x_s2,'linear','extrap');%X_s(2:end)
g=g-mean(s2);

gg = spline(x_s,g,xt);

% dxx = 13100/410;
% xxx = 0:dxx:410-dxx;

figure, plot(gg), hold on, plot(s2-mean(s2),'r')
%%

% exp_name = 'LC1_dt25ms_1';
% path1 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv1\';
% path2 = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\piv2\';
% path = '\\beo\data\ExpLC1_dt25ms_1_WP\RawImages\pivsurf\';
% image_pair_number = str2double(imnum);
% num_of_digits = 4;
% image_letter = 'b';
% % u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
% % fusedIm = fuseImagesLC(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter ,tform);
% 
% % xpiv = 0:pivRes:3784*pivRes;
% % ypiv = 0:pivRes:2048*pivRes;
% % figure, imagesc(xpiv,ypiv,fusedIm.fused_im), colormap(bone), caxis([0 400]), hold on, plot(s2_pivSample-0.3,'r')

% fused_im = 
% u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
% imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% % s3 = imSurfLC.z_s * pivsurfRes;
% imgPivsurf =fliplr(imSurfLC.pivSurf_resized);
% sPS = imSurfLC.z_s(end:-1:1).*pivRes/2;
% % figure,imagesc(imgPivsurf), colormap(bone), hold on, plot(imSurfLC.z_s(end:-1:1), 'r')
% 
% sLfv = s2;
% % figure, plot(s2_pivSample-mean(s2_pivSample)), hold on, plot(s3(end:-1:1)-mean(s3),'r')
% % figure, plot(sLfv-mean(sLfv))
% hold on, plot(sPS-mean(sPS),'k')
% 
% sPSStart = 4900;
% figure, plot(sLfv-mean(sLfv)), hold on, plot(sPSStart:sPSStart+3784,sPS-mean(sPS),'r')
% 
% figure, imagesc(imSurfLC.pivSurf_resized), colormap(bone)






toc













