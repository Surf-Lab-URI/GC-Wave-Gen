clear all
close all

load('\\beo\data\ExpLC4_dt7ms_1\RawImages\Lfv\ExpLC4_dt7ms_1_Lfv_1547.mat')

figure, imagesc(imgLfv), colormap(bone), caxis([0 1000])

%% correct image for lens distortion
load('\\beo\mFiles\lfv\Calib_Results.mat')
I = imgLfv;

KK = [fc(1) alpha_c*fc(1) cc(1) ; 0 fc(2) cc(2) ; 0 0 1];
    
%%% Compute the new KK matrix to fit as much data in the image (in order to
%%% accomodate large distortions:
r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

[I2] = rect(I,eye(3),fc,cc,kc,alpha_c,KK_new);

figure, imagesc(I2), colormap(bone), caxis([0 500])
 
 %% correct image for camera angle
 
U = [ 638 886
     1383 981
     1393 1262
      627 1261];
X = [ 627 886
     1393 981
     1393 1262
      627 1261];
% U = [U(:,2) U(:,1)];
% X = [X(:,2) X(:,1)];
tform = maketform('projective',U,X);
 
% imgPiv1_t=imtransform(imgPiv1,tform,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1); 
imgLfv_t=imtransform(I2,tform,'FillValues',2000);
figure, imagesc(imgLfv_t), colormap(bone), caxis([0 1000])

imgLfv_fin = imgLfv_t(101:2000,173:2092);

lfvSurf = findSurfaceLfv(imgLfv_fin);

lfvSurfI2 = findSurfaceLfv(I2);
figure, plot(lfvSurfI2.z_s(85:2000)), hold on, plot(lfvSurf.z_s-35,'r')
figure, imagesc(imgLfv_fin), colormap(bone), caxis([0 1000]), hold on, plot(lfvSurf.z_s,'r')