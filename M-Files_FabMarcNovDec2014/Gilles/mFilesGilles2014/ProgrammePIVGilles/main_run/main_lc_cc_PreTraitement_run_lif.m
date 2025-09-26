tic,
clear all
close all

%% Parametres
exp_name = 'LC1_4';
ResPivReal = 6.4101e-05; % m/pixel
ResPiv=1;
hResLfv = 4.5663e-4/ResPivReal;
vResLfv = 4.4956e-4/ResPivReal;
num_of_digits = 4;
IntrWndw = [128 64 32 16 8];
GrdSpc = [64 32 16 8 4];
n = 2048; m = 2048;
path = ['E:\data\CamFinale\PivLif\PivLif4\Movie7_Scene7'];
save_path = ['E:\ComputedVelocities\PivLif4\'];
start_image_pair_number = 365;
end_image_pair_number = 365;
% %  Remarque : dt=12.5ms
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]';
y=y-16;
tform = maketform('affine',[x y],[u v]);
epsilon=0.1;
thresh=3;
b=4;
clear u v x y
%%
for image_pair_number = start_image_pair_number:end_image_pair_number
    image_pair_number
    %% Images aux temps 'a'
    image_letter='a';
    %% Chargement de Piv
    rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM1 = matFrame.img;
    IMM1 = medfiltEpsThreshParallele(IM1,epsilon,thresh,b);
    IMM1=inpaint_nans(IMM1,b);
    IM1 = IM1-IMM1;
    %% Chargement et correction de Pivsurf
   
    
    %% Images aux temps 'b'
    image_letter='b';
    %% Chargement de Piv
    rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM2 = matFrame.img;
    IMM2 = medfiltEpsThreshParallele(IM2,epsilon,thresh,b);
    IMM2=inpaint_nans(IMM2,b);
    IM2 = IM2-IMM2;
    %% Chargement et correction de Pivsurf
mask1=ones(size(IM1));
  mask2=ones(size(IM1));  
    
    compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask2, IntrWndw, GrdSpc);
    compVel1 = computeVelocities_Fab(IM1, IM2, mask1, mask2, IntrWndw, GrdSpc);
    
    
    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_pretraitement'];
    outfile = [save_path filename];
    save(outfile, 'compVel','compVel1','SUMASK');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
end
toc