tic,
clear all
close all


%% Parameter
exp_name = 'LC3_dt25ms_1';
deltaT = 7d-3; %s
path = '\\beo\data\';
num_of_digits = 4;
piv_res = 40d-6; %m/pix
vec_res = 4 * piv_res;
pivsurf_res = 102d-6; %m/pix
piv_delta_t = 1/7.2; %sec

pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
iws = 16; % initial widow size
vResLfv =  1.3543e-04/pivResReal; %pixel de piv/pixel de lfv;  %lfv resolution
hResLfv = 1.2912e-04/pivResReal; %pixel de piv/pixel de lfv
lim = 4700; % position of PIVimage in LFV

IntrWndw = [128 64 32 16]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [ 64 32 16 8];

% save_path = '';
% start_image_pair_number = 1150;
% end_image_pair_number = 1150;

u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
u = [ 760 1234 1015]'; v = [ 1131 1128 1093]'; x = [237 3515 1991]'; y = [319 316 74]'; tformLfv = maketform('affine',[u v],[x y]);

clear u v x y
%%
% for image_pair_number = start_image_pair_number:end_image_pair_number
image_pair_number = 790
%% Images aux temps 'a'
image_letter='a';
%% Correction et fusion des images de PIV
piv1_struc = load([path 'Exp' exp_name '\RawImages\Piv1\' 'Exp' exp_name '_Piv1_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
imgPiv1 = fliplr(piv1_struc.imgPiv1);
piv2_struc = load([path 'Exp' exp_name '\RawImages\Piv2\' 'Exp' exp_name '_Piv2_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
imgPiv2 = fliplr(piv2_struc.imgPiv2);
IM1= fuseImages_lc(imgPiv1,imgPiv2,tform);
clear piv1_struc piv2_struc imgPiv1 imgPiv2
%% Chargement et correction de Pivsurf, detection de la surface, creation d'un mask
pivsurf_struc = load([path 'Exp' exp_name '\RawImages\Pivsurf\' 'Exp' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);

ims = pivsurf_struc.imgPivsurf;
ims_L = ims(1:end,1:floor(size(ims,2)/2));
ims_R = ims(1:end,floor(size(ims,2)/2)+1:end);
offs = -9;
ims_fin = [ims_L*1.07 + offs ims_R];
%     ims_L_w = ims_L(1200:end,:)+offs;
%     ims_R_w = ims_R(1200:end,:);
%
%
%     figure, imagesc(ims_L_w), colormap(bone)
%     [nL,xL] = hist(ims_L_w(:),50);
%     [nR,xR] = hist(ims_R_w(:),50);
%     figure, plot(xL*1.07,nL), hold on, plot(xR,nR,'r')
%     figure, imagesc([ims_L*1.07 + offs ims_R]), colormap(bone)
%     offs = -5;
%     ims_L_a = ims_L(1:600,:)+offs;
%     ims_R_a = ims_R(1:600,:);
%     figure, imagesc([ims_L_a*1.07 + offs ims_R_a]), colormap(bone)

imgPivsurf = fliplr(ims_fin);
%     imgPivsurf = fliplr(pivsurf_struc.imgPivsurf);
imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
imgPivsurf_t_cut=imgPivsurf_t(:,:);
surf1=findSurface_lc(imgPivsurf_t_cut);

mask1 = nan(size(imgPivsurf_t_cut));
for i=1:size(imgPivsurf_t_cut,2)
    mask1(round(surf1.z_s(i)):end,i) = 1;
end
%% Chargement et correction de Lfv, detection de la surface, calcul des orbitales
lfv_struc = load([path 'Exp' exp_name '\RawImages\Lfv\' 'Exp' exp_name '_Lfv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
imgLfv = lfv_struc.imgLfv;  %% le fliplr est il necesaire ??
imgLfv_LD = correctLfvLensDist_lc(imgLfv);
imgLfv_t=imtransform(imgLfv_LD,tformLfv);
imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
surfLfv = findSurface_lc(imgLfv_t_cut);

psSWL =  (size(surf1.img,1) - 575 + 1) * pivResReal;
XX = 0:pivResReal:size(IM1,2)*pivResReal-pivResReal;
YYL = size(IM1,1)*pivResReal-pivResReal;
YY = - psSWL : pivResReal : YYL - psSWL;
% YY = fliplr(YY);

%% plot piv fused
figure, imagesc(XX,YY,rot90(IM1,2)), colormap(bone), caxis([200 400])
hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

save_path = 'F:\figures\';
eval(['export_fig ' save_path 'Exp' exp_name '_pivfuse_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])


%% plot pivsurf

figure, imagesc(XX,YY,rot90(surf1.img,2)), colormap(bone)
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

% save_path = 'F:\figures\';
eval(['export_fig ' save_path 'Exp' exp_name '_pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])


%% plot lfv

lfvIm = surfLfv.img;
hResLfv = 1.2912e-04; %m/pix
vResLfv =  1.3543e-04;
% lfvWidth = hResLfv * size(surfLfv.img,2);
lfvSWL = (size(lfvIm,1) - 2295 + 1) * vResLfv; %pixel
% psSWL =  (size(surf1.img,1) - 575 + 1) * pivResReal;
XX = 0:hResLfv:size(lfvIm,2)*hResLfv-hResLfv;
YYL = size(lfvIm,1)*vResLfv-vResLfv;
YY = - lfvSWL : vResLfv : YYL - lfvSWL;

figure, imagesc(XX,YY,flipud(surfLfv.img)), colormap(bone)
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

% save_path = 'F:\figures\';
eval(['export_fig ' save_path 'Exp' exp_name '_lfv_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])