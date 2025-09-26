clear all
close all

%% Ouverture des fichiers
% Image calib
matFramePivsurf_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\pivsurf\Movie11_Scene8_pivsurf_0.raw',2048,2048);
matFrameLfv_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\lfv\Movie11_Scene8_lfv_0.raw',2048,2048);

imgPivsurf=flipud(fliplr(matFramePivsurf_Sc8.img));
imgLfv=flipud(matFrameLfv_Sc8.img);


figure, imagesc(imgLfv), colormap(bone), caxis([])
