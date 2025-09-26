tic,
clear all
close all


%% Parameter
exp_name = 'LC3dt10';
%%
num_of_digits = 4;

IntrWndw = [128 64 32 16 8]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8 4]; 

save_path = ['E:\ComputedVelocities\20140730\LC3dt10\'];
start_image_pair_number = 355;
end_image_pair_number = 355;


u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1088 1091 570]'; 
tform = maketform('affine',[x y],[u v]);


clear u v x y
%%
 for image_pair_number = start_image_pair_number:end_image_pair_number
% image_pair_number = 350

    %% Images aux temps 'a'
    image_letter='a';
     %% Chargement de Piv
     n = 2048;
     m = 2048;
     rawFrame = ['E:\data\20140730\LC3dt10\Movie3_Scene2_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
     matFrame = saveToMatSingleFrame(rawFrame,n,m);
     IM1 = matFrame.img;

     %%Chargement et correction de Pivsurf
     rawFrame = ['E:\data\20140730\LC3dt10\Movie3_Scene2_pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
     matFrame = saveToMatSingleFrame(rawFrame,n,m);
     Pivsurfcc = matFrame.img;
     
     imgPivsurfcc=imtransform(Pivsurfcc,tform,'Xdata',[1 2048],'Ydata',[1 2048], 'FillValues',-1);

     %% Creation du mask
    surf1=findSurface_lc(imgPivsurfcc);
    z_s1 = surf1.z_s;
    mask1 = nan(size(imgPivsurfcc));
    for i=1:size(imgPivsurfcc,2)
        mask1(round(z_s1(i))+10:end,i) = 1; % decalage de 10
    end
    

    
    
    %% Images aux temps 'b'
    image_letter='b';
    %% Correction pivcc et pivsurfcc
    
    rawFrame = ['E:\data\20140730\LC3dt10\Movie3_Scene2_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM2 = matFrame.img;
    
    rawFrame = ['E:\data\20140730\LC3dt10\Movie3_Scene2_pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    Pivsurfcc2 = matFrame.img;
    
    imgPivsurfcc2=imtransform(Pivsurfcc2,tform,'Xdata',[1 2048],'Ydata',[1 2048], 'FillValues',-1);
    
    surf2=findSurface_lc(imgPivsurfcc2);
    z_s2 = surf2.z_s;
    mask2 = nan(size(imgPivsurfcc2));
    for i=1:size(imgPivsurfcc2,2)
        mask2(round(z_s2(i))+10:end,i) = 1; % decalage de 10
    end

    %% piv calc
    compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask2, IntrWndw, GrdSpc);
    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_3232'];
    outfile = [save_path filename];
    save(outfile, 'compVel');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
 end
toc