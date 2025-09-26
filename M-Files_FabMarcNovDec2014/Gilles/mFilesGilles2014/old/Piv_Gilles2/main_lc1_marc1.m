tic,
clear all
close all


%% Parameters
exp_name = 'LC2_dt7ms_1';
deltaT = 7d-3; %s
path = '\\beo\data\';
num_of_digits = 4;

pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
iws = 8; % initial window size

IntrWndw = [512 256 128 64 32 16];
GrdSpc = [256 128 64 32 16 8];

save_path = '';
start_image_pair_number = 1150;
end_image_pair_number = 1150;

u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
clear u v x y
%%
% for image_pair_number = start_image_pair_number:end_image_pair_number
image_pair_number = 1400
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
    tic
    %% Calcul des orbitales sur surf1
%     sPSa = surf1.z_s(end:-1:1).*pivRes;
%     sPS = size(surf1.img,1)*pivRes - sPSa;
%     Orbitals = calcOrbitals_lc(sPS,pivRes,iws,deltaT,pivResReal);
%     Orbitals_interp = interpOrbitals_lc(sPS,pivRes,Orbitals.u,Orbitals.w,iws);
%     toc
%     %% Calcul de la PIV
% %     tic
%     delxOrb = Orbitals_interp.u;
%     delzOrb = Orbitals_interp.w;
%     tic
%     c2 = PIV_FAB6_PivsurfOrb1 (IM1,IM2, mask1, mask2, delxOrb, delzOrb, IntrWndw, GrdSpc, iws);
%     toc
%     compVel_orb1 = PIV_FAB6_PivsurfOrb1 (IM1,IM2, mask1, mask2, delxOrb, delzOrb, IntrWndw, GrdSpc, iws);
%     toc
% end

%%
mask12 = mask1 + mask2;
% compVel =  computeVelocities_LC1(fliplr(IM1), fliplr(IM2), fliplr(mask12));
compVel_1 = PIV_FAB6_PivsurfOrb1_marc1 (IM1,IM2, mask12, mask12,IntrWndw, GrdSpc);

close all
% figure, imagesc(compVel.delta_x.* compVel.mask), caxis([-20 37])
figure, imagesc(compVel_1.delx_int.*compVel_1.MASK)%, caxis([-20 37])
figure, imagesc(compVel_1.dely_int.*compVel_1.MASK)%, caxis([-20 37])

% figure, imagesc(compVel.cor_mtx), caxis([0 1])
figure, imagesc(compVel_1.dcor)%, caxis([0 1])

% dcor
numel(find(compVel_1.dcor<0.5))/numel(compVel_1.dcor)

% %% testing results
dcor = compVel_1.dcor;
MASK = compVel_1.MASK;
dcor1 = 50 * ones(size(dcor)); % matrix of 50s
dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
figure, imagesc(dcor1), caxis([0 1]), colorbar
numel(find(dcor<0.5))/numel(find(dcor1<40))
% 
% mask111 = mask11;
% mask222 = mask22;
% mask111(isnan(mask11)) = 50;
% mask222(isnan(mask22)) = 50;
% mask_diff = mask222-mask111;
% figure, imagesc(mask111)
% figure, imagesc(mask_diff), colorbar

toc