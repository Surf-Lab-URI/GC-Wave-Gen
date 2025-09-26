clear all
close all

hRes = 2.6336d-4; %m/pixel
vRes = 2.7152d-4; %m/pixel
pivRes = 30d-6; %m/pixel
lfvToPivScale = hRes/pivRes;

rawFrame = '\\beo\data\rawPIV\Movie11_Scene8\lfv\Movie11_Scene8_lfv_0.raw';
n = 2048;
m = 2048;
calibMatFrame = saveToMatSingleFrame(rawFrame,n,m);
imgLfv = calibMatFrame.img;
figure, imagesc(imgLfv), colormap(bone), caxis([0 1000])

% load('\\beo\data\ExpLC4_dt7ms_1\RawImages\Lfv\ExpLC4_dt7ms_1_Lfv_1547.mat')

imgLfv_LD = correctLfvLensDist(imgLfv);
imgLfv_fin = correctLfvCamAngle(imgLfv_LD);

figure, imagesc(flipud(imgLfv_fin)), colormap(bone), caxis([0 1000])

%%
rawFrame = '\\beo\data\rawPIV\Movie11_Scene8\piv1\Movie11_Scene8_piv1_0.raw';
n = 2048;
m = 2048;
calibMatFrame = saveToMatSingleFrame(rawFrame,n,m);
imgPiv1 = calibMatFrame.img;
figure, imagesc(imgPIv1), colormap(bone), caxis([0 1000])
