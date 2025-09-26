clear all
close all
rawFrame = '\\beo\data\rawPIV\Movie11_Scene8\lfv\Movie11_Scene8_lfv_0.raw';
n = 2048;
m = 2048;
calibMatFrame = saveToMatSingleFrame(rawFrame,n,m);

figure, imagesc(calibMatFrame.img), colormap(bone), caxis([0 1000])

%% correct image for lens distortion
load('\\beo\mFiles\lfv\Calib_Results.mat')
I = calibMatFrame.img;

KK = [fc(1) alpha_c*fc(1) cc(1) ; 0 fc(2) cc(2) ; 0 0 1];
    
    %%% Compute the new KK matrix to fit as much data in the image (in order to
    %%% accomodate large distortions:
    r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
    dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
    fc_new = dist_amount * fc;
    
    KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];
    
    [I2] = rect(I,eye(3),fc,cc,kc,alpha_c,KK_new);
 %% quick lfv resolution calculation
%  figure, imagesc(I2), colormap(bone), caxis([0 1000])
%  hRes = ginput(2);
%  hRes = (9*2.54d-2)/(hRes(2,1) - hRes(1,1)) %m/pixel
%  vRes = ginput(2);
%  vRes = (3.5*2.54d-2)/(vRes(2,2) - vRes(1,2)) %m/pixel
%  %%
%  hRes = 2.6941d-4; %m/pixel
%  vRes = 2.7098d-4; %m/pixel
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
imgLfv_t = imtransform(I2,tform,'FillValues',2000); 
figure, imagesc(imgLfv_t), colormap(bone), caxis([0 1000])

hRes = ginput(2);
hRes = (9*2.54d-2)/(hRes(2,1) - hRes(1,1)) %m/pixel
hRes = 2.6336d-4; %m/pixel

vRes = ginput(2);
vRes = (4*2.54d-2)/(vRes(2,2) - vRes(1,2)) %m/pixel
vRes = 2.7152d-4; %m/pixel








