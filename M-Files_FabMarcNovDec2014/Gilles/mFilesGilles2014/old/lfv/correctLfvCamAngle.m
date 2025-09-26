function imgLfv_CA = correctLfvCamAngle(imgLfv_LD)

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
imgLfv_t=imtransform(imgLfv_LD,tform,'FillValues',2000);
% figure, imagesc(imgLfv_t), colormap(bone), caxis([0 1000])

imgLfv_CA = imgLfv_t(101:2000,173:2092);

% lfvSurf = findSurfaceLfv(imgLfv_fin);

% lfvSurfI2 = findSurfaceLfv(I2);
% figure, plot(lfvSurfI2.z_s(85:2000)), hold on, plot(lfvSurf.z_s-35,'r')
% figure, imagesc(imgLfv_fin), colormap(bone), caxis([0 1000]), hold on, plot(lfvSurf.z_s,'r')