clear all
close all

matFrameLfv_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\lfv\Movie11_Scene8_lfv_0.raw',2048,2048);
imgLfv = matFrameLfv_Sc8.img;

imgLfv_t = correctLfv1(imgLfv);

imgLfv_t_cut = imgLfv_t(21:4077, 158:4216);

figure, imagesc(imgLfv_t_cut), colormap(bone), caxis([0 1000])

vCal = ginput(2);
vRes = (3*2.54e-2)/(vCal(2,2) - vCal(1,2));

vRes =  1.3543e-04; %m/pix;

hCal = ginput(2);
hRes = (9*2.54e-2)/(hCal(2,1) - hCal(1,1));
hRes = 1.2912e-04; %m/pix;