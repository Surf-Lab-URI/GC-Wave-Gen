tic,
clear all
close all

%% Parametres
exp_name = 'LC2_2';
ResPivReal = 6.4101e-05; % m/pixel
ResPiv=1;
hResLfv = 4.5663e-4/ResPivReal; 
vResLfv = 4.4956e-4/ResPivReal; 
num_of_digits = 4;
IntrWndw = [128 64 32 16 8];
GrdSpc = [64 32 16 8 4];
n = 2048; m = 2048;
path = ['E:\data\20140804\LC2_2\Movie4_Scene6'];
save_path = ['E:\ComputedVelocities\20140804\LC2_2\'];
start_image_pair_number = 340;
end_image_pair_number = 420;
% %  Remarque : dt=12.5ms
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]'; 
y=y-16;
tform = maketform('affine',[x y],[u v]);
clear u v x y
%%
 for image_pair_number = start_image_pair_number:end_image_pair_number
% image_pair_number = 365
image_pair_number
    %% Images aux temps 'a'
    image_letter='a';
     %% Chargement de Piv
     rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
     matFrame = saveToMatSingleFrame(rawFrame,n,m);
     IM1 = matFrame.img;
     %% Chargement et correction de Pivsurf
     rawFrame = [path '_pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
     matFrame = saveToMatSingleFrame(rawFrame,n,m);
     Pivsurfcc = matFrame.img;
     imgPivsurfcc=imtransform(Pivsurfcc,tform,'Xdata',[1 2048],'Ydata',[1 2048], 'FillValues',-1);
     M =calcM(imgPivsurfcc);
     clear image_letter rawFrame matFrame Pivsurfcc
     %% Creation du mask
    surf1=findSurface_lc(M);
    mask1 = nan(size(M));
    for i=1:size(M,2)
        mask1(round(surf1.z_s_f(i))+0:end,i) = 1; % decalage de 0
    end
   
    %% Images aux temps 'b'
    image_letter='b';
    %% Chargement de Piv
    rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM2 = matFrame.img;
    %% Chargement et correction de Pivsurf
    rawFrame = [path '_pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    Pivsurfcc2 = matFrame.img;
    imgPivsurfcc2=imtransform(Pivsurfcc2,tform,'Xdata',[1 2048],'Ydata',[1 2048], 'FillValues',-1);
    M2 =calcM(imgPivsurfcc2);
    
    clear image_letter rawFrame matFrame Pivsurfcc2
    %% Creation du mask
    surf2=findSurface_lc(M2);
    mask2 = nan(size(M2));
    for i=1:size(M2,2)
        mask2(round(surf2.z_s_f(i))+0:end,i) = 1; % decalage de 0
    end

    %% Creation de SU 
    rawFrame = [path '_lfv_' sprintf(['%0' num2str(num_of_digits-1) 'd'], image_pair_number) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    Lfv= rot90(matFrame.img,2);
    Lfv_LD = correctLfvLensDist_lc(Lfv);
    Lfv_cut=Lfv_LD(401:1600,401:1600); % Seul une partie de l'image est expoitable
    surfLfv=findSurface_lc(Lfv_cut); 
    
    %% Calcul de la PIV
    SU1 = SUcc(surfLfv, surf1, ResPiv, hResLfv, vResLfv, GrdSpc(end),1);
    SU2 = SUcc(surfLfv, surf2, ResPiv, hResLfv, vResLfv, GrdSpc(end),1);
    dSU = SU2-SU1;  %% Je serais tenter de rajouter un gros smoothn bien costaud N.B. Le faire sur les surf1 et surf2, moins couteux
%     Interpoler dSU sur surf1 et s'en servir comme initialisation pour
%     compVel
    compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask2, IntrWndw, GrdSpc);
    compVel1 = computeVelocities_Fab(IM1, IM2, mask1, mask2, IntrWndw, GrdSpc);
   
    SurfMASK.img=compVel.MASK;
    SurfMASK.z_s_f = nansum(compVel.MASK);  % On recupere la surface
    SurfMASK.z_s_f = size(compVel.MASK,1) - SurfMASK.z_s_f; % On la remet dans "le bon sens"
    SUMASK = SUcc(surfLfv, SurfMASK, ResPiv, hResLfv, vResLfv, 1,1);
    dSUInterp = interpTransfoNew(dSU,SUMASK);
    delycor=compVel.dely_ints+dSUInterp;
    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel','compVel1','SUMASK');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
 end
toc