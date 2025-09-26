%% Calcul les orbitales hybrides et les enregistre a la place des precedentes


tic,
clear all
close all


%% Parameter
exp_name = 'LC4_dt25ms_1';
deltaT = 25d-3; %s

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
lim_left_big=697; % position of PIVimage in PIVsurf (when it's not resized to PIV - i.e. pivsurfbig)

IntrWndw = [128 64 32 16]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8];

save_path = ['\\beo\data\Exp' exp_name '\ComputedVelocities\'];
start_image_number = 299;
end_image_number =299;

u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
u = [ 760 1234 1015]'; v = [ 1131 1128 1093]'; x = [237 3515 1991]'; y = [319 316 74]'; tformLfv = maketform('affine',[u v],[x y]);


for image_number = start_image_number:end_image_number
    image_letter='a';
    
    pivsurf_struc = load([path 'Exp' exp_name '\RawImages\Pivsurf\' 'Exp' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number) '_' image_letter]);
    imgPivsurf = fliplr(pivsurf_struc.imgPivsurf);
    imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
    imgPivsurf_big=imtransform(imgPivsurf,tformSurf,'XYScale',1);
    
    imgPivsurf_t_cut=imgPivsurf_t(:,:);
    surf1=findSurface_lc(imgPivsurf_t_cut);
    lim_right=length(surf1.img);
    surfbig=findSurface_lc(imgPivsurf_big);
    mask1 = nan(size(imgPivsurf_t_cut));
    for i=1:size(imgPivsurf_t_cut,2)
        mask1(round(surf1.z_s(i)):end,i) = 1;
    end
    
    lfv_struc = load([path 'Exp' exp_name '\RawImages\Lfv\' 'Exp' exp_name '_Lfv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number)]);
    imgLfv = lfv_struc.imgLfv;
    imgLfv_LD = correctLfvLensDist_lc(imgLfv);
    imgLfv_t=imtransform(imgLfv_LD,tformLfv);
    imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);
    surfLfv = findSurface_lc(imgLfv_t_cut);
    
    Orbitals = calcOrbitals_hybrid(surfLfv,surfbig,vResLfv,hResLfv,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right);
    Surface = calcSurface_lc(surfbig,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_right);
    Orbitals_interp = interpOrbitals_lc_lfv2(surf1,pivRes,Orbitals.u,Orbitals.w,Surface.SU,IntrWndw,GrdSpc);
    
    load([save_path 'Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number)]);
    
    compVel.delxOrb= Orbitals_interp.u;
    compVel.delyOrb= Orbitals_interp.w;

    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel');
    disp(['pair ' num2str(image_number) ' done.']);

    
end