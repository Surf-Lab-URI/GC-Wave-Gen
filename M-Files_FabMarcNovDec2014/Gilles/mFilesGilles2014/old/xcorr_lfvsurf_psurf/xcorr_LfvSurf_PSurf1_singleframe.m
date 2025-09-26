tic,
clear all
close all


%% Parameters
exp_name = 'LC3_dt7ms_1';
deltaT = 7d-3; %s
path = '\\beo\data\';
num_of_digits = 4;

pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
iws = 8; % initial window size

IntrWndw = [512 256 128 64 32 16];
GrdSpc = [256 128 64 32 16 8];

save_path = '';
start_image_pair_number = 1000;
end_image_pair_number = 2000;
lag_vec = nan(end_image_pair_number - start_image_pair_number + 1,1);

u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
u = [ 760 1234 1015]'; v = [ 1131 1128 1093]'; x = [237 3515 1991]'; y = [319 316 74]'; tformLfv = maketform('affine',[u v],[x y]);
clear u v x y
%%
% for image_pair_number = start_image_pair_number:end_image_pair_number
     image_pair_number = 1034
    disp(['pair ' num2str(image_pair_number)]);
   
    %% Images aux temps 'a'
    image_letter='a';
    
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
    
    
    
    lfv_struc = load([path 'Exp' exp_name '\RawImages\Lfv\' 'Exp' exp_name '_Lfv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
    imgLfv = lfv_struc.imgLfv;  %% le fliplr est il necesaire ??
    imgLfv_LD = correctLfvLensDist_lc(imgLfv);
    imgLfv_t=imtransform(imgLfv_LD,tformLfv);
    imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
    surfLfv = findSurface_lc(imgLfv_t_cut);
    
%     figure, imagesc(surf1.img), colormap(bone), hold on, plot(surf1.z_s,'r')
%     figure, imagesc(surfLfv.img), colormap(bone), hold on, plot(surfLfv.z_s,'r')
%     
    
    %% interpolate Lfv surface to piv grid
    
    piv_res = 40d-6; %m/pix
    psSWL =  size(surf1.img,1) - 575 + 1; %pixel
    lfvSWL = size(surfLfv.img,1) - 2295 + 1; %pixel
    hResLfv = 1.2912e-04; %m/pix
    vResLfv =  1.3543e-04;
    lfvWidth = hResLfv * size(surfLfv.img,2);
    
    z_s_lfv = (size(surfLfv.img,1) - surfLfv.z_s + 1 - lfvSWL) * vResLfv;
    z_s_ps = (size(surf1.img,1) - surf1.z_s + 1 - psSWL) * piv_res;
    
    % z_s_lfv = z_s_lfv(end:-1:1);
    z_s_ps = z_s_ps(end:-1:1);
    
    %%
    
    xi = 0 : hResLfv : lfvWidth - hResLfv;
    xt = 0 : piv_res : lfvWidth - piv_res;
    z_s_lfv_interp = spline(xi,z_s_lfv,xt);
    figure, plot(z_s_lfv_interp), hold on, plot(4700:4700-1+length(z_s_ps),z_s_ps,'r')
    %% plot surface from pivsurf (piv_fused)
    
    
    % figure, plot(z_s_ps)
    
    %% find position of piv_fused within lfv
    [c,lags] = xcorr(z_s_lfv_interp, z_s_ps);
    c(lags<3700) = 0;
    c(lags>5700) = 0;
    figure, plot(lags,c);
    [M,Ind] = max(c);
    lag_vec(image_pair_number) = lags(Ind);
    figure, plot(z_s_lfv_interp), hold on, plot(lags(Ind):lags(Ind)-1+length(z_s_ps),z_s_ps,'r')
%     size(surfLfv.img)
% end
%%
% save('\\beo\mFiles\xcorr_lfvsurf_psurf\lag_vec.mat','lag_vec')
toc