tic,
clear all
close all


%% Parameter
exp_name = 'LC2_dt7ms_1';
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

IntrWndw = [512 256 128 64 32 16]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [512 256 64 32 16 8]; 

save_path = '';
start_image_pair_number = 1150;
end_image_pair_number = 1150;

u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
u = [ 760 1234 1015]'; v = [ 1131 1128 1093]'; x = [237 3515 1991]'; y = [319 316 74]'; tformLfv = maketform('affine',[u v],[x y]);

clear u v x y
%%
% for image_pair_number = start_image_pair_number:end_image_pair_number
image_pair_number = 1150
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
    imgPivsurf = fliplr(pivsurf_struc.imgPivsurf);
    imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
    imgPivsurf_t_cut=imgPivsurf_t(:,:);
    surf1=findSurface_lc(imgPivsurf_t_cut);
    mask1 = nan(size(imgPivsurf_t_cut));
    for i=1:size(imgPivsurf_t_cut,2)
        mask1(round(surf1.z_s(i)):end,i) = 1;
    end
    clear pivsurf_struc imgPivsurf imgPivsurf_t_cut
    %% Images aux temps 'b'
    image_letter='b';
    %% Correction et fusion des images de PIV
    piv1_struc = load([path 'Exp' exp_name '\RawImages\Piv1\' 'Exp' exp_name '_Piv1_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    imgPiv1 = fliplr(piv1_struc.imgPiv1);
    piv2_struc = load([path 'Exp' exp_name '\RawImages\Piv2\' 'Exp' exp_name '_Piv2_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    imgPiv2 = fliplr(piv2_struc.imgPiv2);
    IM2= fuseImages_lc(imgPiv1,imgPiv2,tform);
    clear piv1_struc piv2_struc imgPiv1 imgPiv2
    %% Chargement et correction de Pivsurf, detection de la surface, creation d'un mask
    pivsurf_struc = load([path 'Exp' exp_name '\RawImages\Pivsurf\' 'Exp' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    imgPivsurf = fliplr(pivsurf_struc.imgPivsurf);
    imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
    imgPivsurf_t_cut=imgPivsurf_t(:,:);
    surf2=findSurface_lc(imgPivsurf_t_cut);
    mask2 = nan(size(imgPivsurf_t_cut));
    for i=1:size(imgPivsurf_t_cut,2)
        mask2(round(surf2.z_s(i)):end,i) = 1;
    end
    clear pivsurf_struc imgPivsurf imgPivsurf_t_cut image_letter
    %% Chargement et correction de Lfv, detection de la surface, calcul des orbitales
    lfv_struc = load([path 'Exp' exp_name '\RawImages\Lfv\' 'Exp' exp_name '_Lfv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
    imgLfv = lfv_struc.imgLfv;  %% le fliplr est il necesaire ?? 
    imgLfv_LD = correctLfvLensDist_lc(imgLfv);
    imgLfv_t=imtransform(imgLfv_LD,tformLfv);
    imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
    surfLfv = findSurface_lc(imgLfv_t_cut);
    Orbitals = calcOrbitals_lc_lfv(surfLfv,imgLfv_t_cut,vResLfv,hResLfv,pivRes,iws,deltaT,pivResReal,IntrWndw,GrdSpc,lim,length(surf1.img));
    Orbitals_interp = interpOrbitals_lc_lfv(surf1,pivRes,Orbitals.u,Orbitals.w,IntrWndw,GrdSpc);
    %% Calcul de la PIV
    delxOrb=Orbitals_interp.u;
    delzOrb=Orbitals_interp.w;
    delxOrb1 = flipud(delxOrb);
    delzOrb1 = flipud(delzOrb);
    IM11 = fliplr(IM1);
    IM22 = fliplr(IM2);
    mask11 = fliplr(mask1);
    mask22 = fliplr(mask2);
    compVel_orb2 = PIV_FAB6_LfvOrb1 (IM11,IM22, mask11, mask22, delxOrb1, delzOrb1, IntrWndw, GrdSpc);
% end

%% testing results

dcor = compVel_orb2.dcor;
MASK = compVel_orb2.MASK;
dcor1 = 50 * ones(size(dcor)); % matrix of 50s
dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
figure, imagesc(dcor1), caxis([0 1]), colorbar
numel(find(dcor<0.5))/numel(find(dcor1<40))

%%
Turbx=compVel_orb2.delx_int(1:255,1:472)-compVel_orb2.delxOrb(1:255,1:472);
Turby=compVel_orb2.dely_int(1:255,1:472)-compVel_orb2.delyOrb(1:255,1:472);
figure,imagesc(Turbx)
figure,imagesc(Turby)

% figure, quiver(1:10:472,1:5:255,compVel_orb2.delx_int(1:5:255,1:10:472).* compVel_orb2.MASK(1:5:255,1:10:472),compVel_orb2.dely_int(1:5:255,1:10:472).* compVel_orb2.MASK(1:5:255,1:10:472))



