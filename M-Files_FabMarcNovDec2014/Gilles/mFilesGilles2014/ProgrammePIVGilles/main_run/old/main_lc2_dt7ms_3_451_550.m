tic,
clear all
close all


%% Parameter
exp_name = 'LC2_dt7ms_3';
deltaT = 7d-3; %s

%%
path = '\\beo\data\';
num_of_digits = 4;
piv_res = 40d-6; %m/pix
vec_res = 4 * piv_res;
pivsurf_res = 102d-6; %m/pix
piv_delta_t = 1/7.2; %sec

pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
vResLfv =  1.3543e-04/pivResReal; %pixel de piv/pixel de lfv;  %lfv resolution
hResLfv = 1.2912e-04/pivResReal; %pixel de piv/pixel de lfv
lim_left = 4700; % position of PIVimage in LFV
lim_left_big = 697; % position of PIVimage in PIVsurf (when it's not resized to PIV - i.e. pivsurfbig)

IntrWndw = [256 128 64 32 16]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [128 64 32 16 8]; 

save_path = ['\\beo\data\Exp' exp_name '\ComputedVelocities\'];
start_image_pair_number = 451;
end_image_pair_number = 550;

u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
u = [ 760 1234 1015]'; v = [ 1131 1128 1093]'; x = [237 3515 1991]'; y = [319 316 74]'; tformLfv = maketform('affine',[u v],[x y]);

clear u v x y
%%
 for image_pair_number = start_image_pair_number:end_image_pair_number
% image_pair_number = 612
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
    imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);%pixel PIV resolution cut to PIV size to make mask
    imgPivsurf_big=imtransform(imgPivsurf,tformSurf,'XYScale',1); %pixel PIV resolution not cut to PIV image to avoid edge effects?

    imgPivsurf_t_cut=imgPivsurf_t(:,:);
    surf1=findSurface_lc(imgPivsurf_t_cut);
    lim_right=length(surf1.img);
    surfbig=findSurface_lc(imgPivsurf_big);
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
    imgLfv_t=imtransform(imgLfv_LD,tformLfv); % took out XSCALE to avoid enormous image
    imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
    surfLfv = findSurface_lc(imgLfv_t_cut);
    Orbitals = calcOrbitals_hybrid(surfLfv,surfbig,vResLfv,hResLfv,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right);
%     Orbitals = calcOrbitals_lfv2(surfLfv,surfbig,vResLfv,hResLfv,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right);
%     Surface = calcSurface_lc(surfbig,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_right); %deja calcule avec calcObrbitals-hybrid

%     Ajouter la possibilite de lire le courant moyen en fonction de la profondeur
%     calcule au temps n-1 a Orbitals.u pour une meilleure
%     initialisation
%     for i=1,size(Orbitals.u,1)
%         Orbitals.u(i,:)=Orbitals.u(i,:)+Courant(i);
%     end
    
    Orbitals_interp = interpOrbitals_lc_lfv2(surf1,pivRes,Orbitals.u,Orbitals.w,Orbitals.SU,IntrWndw,GrdSpc);
    %% Calcul de la PIV
    delxOrb=Orbitals_interp.u;
    delzOrb=Orbitals_interp.w;
      
    
    IM11 = fliplr(IM1);
    IM22 = fliplr(IM2);
    mask11 = fliplr(mask1);
    mask22 = fliplr(mask2);
    compVel = PIV_FAB6_LfvOrb2 (IM11,IM22, mask11, mask22, delxOrb, delzOrb, IntrWndw, GrdSpc);
    
%%     Calcul du courant moyen
    imgTransfo255472x = interpTransfo255472(surf1,compVel.delx_ints,Orbitals.SU,GrdSpc);
    imgTransfo255472y = interpTransfo255472(surf1,compVel.dely_ints,Orbitals.SU,GrdSpc);
    for i=1:size(imgTransfo255472x,1)
        Courantx(i)=mean(imgTransfo255472x(i));
        Couranty(i)=mean(imgTransfo255472y(i));
    end
    compVel.courantx=Courantx;
    compVel.couranty=Couranty;
    
%%     Sauvegarde
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    


% 
% dcor = compVel.dcor;
% MASK = compVel.MASK;
% dcor1 = 50 * ones(size(dcor)); % matrix of 50s
% dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
% figure, imagesc(dcor1), caxis([0 1]), colorbar
% numel(find(dcor<0.5))/numel(find(dcor1<40))


 end
% 

toc
%% testing results
% dcor = compVel_orb2.dcor;
% MASK = compVel_orb2.MASK;
% dcor1 = 50 * ones(size(dcor)); % matrix of 50s
% dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
% figure, imagesc(dcor1), caxis([0 1]), colorbar
% 
% mask111 = mask11;
% mask222 = mask22;
% mask111(isnan(mask11)) = 50;
% mask222(isnan(mask22)) = 50;
% mask_diff = mask222-mask111;
% figure, imagesc(mask111)
% figure, imagesc(mask_diff), colorbar

%toc

%%
% figure,imagesc(compVel.delx_int.*compVel.MASK), colorbar, caxis([-25 25])
% figure,imagesc(compVel.dely_int.*compVel.MASK), colorbar, caxis([-25 25])
% figure,imagesc(compVel.dwhatlevel), colorbar
% 
% 
% dcor = compVel.dcor;
% MASK = compVel.MASK;
% dcor1 = 50 * ones(size(dcor)); % matrix of 50s
% dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
% figure, imagesc(dcor1), caxis([0 1]), colorbar
% numel(find(dcor<0.5))/numel(find(dcor1<40))