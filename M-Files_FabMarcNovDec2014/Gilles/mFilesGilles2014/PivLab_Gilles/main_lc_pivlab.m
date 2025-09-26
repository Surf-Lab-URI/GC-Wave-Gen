tic,
clear all
close all

%% Parametres
exp_name = 'LC1_1';
ResPivReal = 6.4101e-05; % m/pixel
ResPiv=1;
hResLfv = 4.5663e-4/ResPivReal; 
vResLfv = 4.4956e-4/ResPivReal; 
num_of_digits = 4;
IntrWndw = [128 64 32 16 8];
GrdSpc = [64 32 16 8 4];
n = 2048; m = 2048;
path = ['E:\data\20140804\LC1_1\Movie4_Scene1'];
save_path = ['E:\ComputedVelocities\20140804\LC1_1\'];
start_image_pair_number = 355;
end_image_pair_number = 380;
% %  Remarque : dt=12.5ms
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]'; 
y=y-16;
tform = maketform('affine',[x y],[u v]);
clear u v x y
%%
%  for image_pair_number = start_image_pair_number:end_image_pair_number
image_pair_number = 368;
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
    mask1 = ones(size(M));
    for i=1:size(M,2)
        mask1(round(surf1.z_s_f(i))+0:end,i) = 0; % decalage de 0
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
    mask2 = ones(size(M2));
    for i=1:size(M2,2)
        mask2(round(surf2.z_s_f(i))+0:end,i) = 0; % decalage de 0
    end

    
    [xtable, ytable, utable, vtable, typevector] = compVel_Pivlab (IM1,IM2,64, 32, 2, mask1, [1 1 2047 2047],4,32,16,8,'linear');
    
    filename = ['Exp' exp_name '_compVel_PivLab_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'xtable','ytable','utable','vtable','typevector');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
%  end
    toc