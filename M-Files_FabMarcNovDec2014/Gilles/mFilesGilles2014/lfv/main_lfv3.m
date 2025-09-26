tic,
clear all
close all

pivRes = 40d-6; %m/pixel
iws = 32; % initial widow size

imnum = '1150';
lim = 4700;

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
z_s = 0 : pivRes * iws : 10d-2; %vertical coordinate up to 20cm with 1mm resolution (to be adjusted to match PIV)
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

%%
exp_name = 'LC2_dt7ms_1';
path1 = ['\\beo\data\Exp' exp_name '\RawImages\Piv1\'];
path2 = ['\\beo\data\Exp' exp_name '\RawImages\Piv2\'];
path = ['\\beo\data\Exp' exp_name '\RawImages\Pivsurf\'];
image_pair_number = str2double(imnum);
num_of_digits = 4;
image_letter = 'a';

u1 = [ 598 1577 990]'; v1 = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u1 v1],[x y]);
imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% s3 = imSurfLC.z_s * pivsurfRes;
imgPivsurf = fliplr(imSurfLC.pivSurf_resized);
sPSa = imSurfLC.z_s(end:-1:1).*pivRes;
%%
image_letter = 'b';
imSurfLC = findSurfaceLC(exp_name, path, image_pair_number, num_of_digits, image_letter,tformSurf);
% s3 = imSurfLC.z_s * pivsurfRes;
imgPivsurf = fliplr(imSurfLC.pivSurf_resized);
sPSb = imSurfLC.z_s(end:-1:1).*pivRes;
% figure,imagesc(imgPivsurf), colormap(bone), hold on, plot(imSurfLC.z_s(end:-1:1), 'r')
% 

imPS = imSurfLC.pivSurf_resized;
sPS = sPSa;
sPS = size(imPS,1)*pivRes - sPS;
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
%     wwi(1:length(xii),col) = spline(xii,ww(1:length(xii),col),xxt);
end

figure, imagesc(uui), caxis([-0.1 0.1])

% figure, imagesc(u), colorbar, caxis([-0.2 0.2])
% figure, imagesc(w), colorbar, caxis([-0.2 0.2])
% beg = 1;
% figure, quiver(beg:10:409, beg:5:157, u(beg:5:157,beg:10:409), w(beg:5:157,beg:10:409))
% 
% %%
% %Surface RECONSTRUCTION at t
% x_s2 = x_s;%0:13100;
% dx2 = dx;
% t = 7d-3;
% f = fft(s2);
% fa = 2*abs(f)/length(f);
% fp = angle(f);
% k = [0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx2;
% omega = sqrt(9.81*k);
% g = 0;
% for i = 1:floor(length(s2)/2)
%     g = g+(fa(i))*cos(k(i)*x_s2(1:end-1)+fp(i)-omega(i)*t)-mean(s2)/floor(length(s2)/2);
% end
% g = interp1(x_s2(1:end-1),g,x_s2,'linear','extrap');%X_s(2:end)
% g = g-mean(s2);
% 
% gg = spline(x_s,g,xt);
% 
% % dxx = 13100/410;
% % xxx = 0:dxx:410-dxx;
% s2m = s2-mean(s2);
% figure, plot(gg), hold on, plot(s2m),'r')


%%


% % figure, plot(sPSa), hold on, plot(sPS,'r')
% 
% sLfv = s2;
% % % figure, plot(s2_pivSample-mean(s2_pivSample)), hold on, plot(s3(end:-1:1)-mean(s3),'r')
% % % figure, plot(sLfv-mean(sLfv))
% % hold on, plot(sPS-mean(sPS),'k')
% % 
% sPSStart = 4700;
% figure, plot(sLfv-mean(sLfv)), hold on, plot(sPSStart:sPSStart+3784,sPSa-mean(sPSa),'r')
% figure, plot(sLfv), hold on, plot(sPSStart:sPSStart+3784,sPSa,'r')
% % 
% % figure, imagesc(imSurfLC.pivSurf_resized), colormap(bone)
% 
% 
% 
% figure, plot(gg(lim:lim+3784)), hold on, plot(sPSb-mean(sPSb),'r')
% figure, plot(s2m(lim:lim+3784)), hold on, plot(sPSa-mean(sPSa),'r')
% % figure, plot(gg(lim:lim+3784)), hold on, plot(s2m(lim:lim+3784),'r')


toc













