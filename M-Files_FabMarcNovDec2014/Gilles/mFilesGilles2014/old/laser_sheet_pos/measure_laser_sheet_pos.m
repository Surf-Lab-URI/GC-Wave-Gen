% measure laser sheet positions cc
close all
clear all
n = 2048;
m = 2048;

% rawFrame = ['F:\data\rawPIV\cc_laser_sheet_pos\Movie13_Scene32_pivsurfcc_0.raw'];
rawFrame = ['F:\data\Calib_cc\Movie13_Scene1_pivsurfcc_0.raw'];
pivsurf_cc = saveToMatSingleFrame(rawFrame,n,m);
pivsurf_cc.img =flipud(fliplr(pivsurf_cc.img));
figure, imagesc(pivsurf_cc.img), colormap(bone), caxis([0 500])

% rawFrame = ['F:\data\rawPIV\cc_laser_sheet_pos\Movie13_Scene32_pivcc_0.raw'];
rawFrame = ['F:\data\Calib_cc\Movie13_Scene1_pivcc_0.raw'];
piv_cc = saveToMatSingleFrame(rawFrame,n,m);
figure, imagesc(piv_cc.img), colormap(bone), caxis([0 1500])

% rawFrame = ['F:\data\rawPIV\cc_laser_sheet_pos\Movie13_Scene30_pivsurf_0.raw'];
% pivsurf = saveToMatSingleFrame(rawFrame,n,m);
% pivsurf.img = flipud(pivsurf.img);
% 
% rawFrame = ['F:\data\rawPIV\cc_laser_sheet_pos\Movie13_Scene30_lfv_0.raw'];
% lfv = saveToMatSingleFrame(rawFrame,n,m);


figure, imagesc(piv_cc.img), colormap(bone)

figure, plot(lfv.img(600,:))
figure, plot(pivsurf.img(600,:))
figure, plot(piv_cc.img(760,:))
figure, plot(pivsurf_cc.img(1070,:))


pscc = pivsurf_cc.img;
% pscc(1070,780) = 100000;
% figure, imagesc(pscc), colormap(bone), caxis([0 50000])

pvcc = piv_cc.img;
% pvcc(760,715) = 100000;


load('F:\data\ExpLC2_dt25ms_cc_1\RawImages\Pivcc\ExpLC2_dt25ms_cc_1_Pivcc_0299_a.mat')
load('F:\data\ExpLC2_dt25ms_cc_1\RawImages\Pivsurfcc\ExpLC2_dt25ms_cc_1_Pivsurfcc_0299_a.mat')

u = [231 1814 204 1855]'; v = [753 753 1913 1913]'; x = [231 1814 231 1814]'; y =[753 753 1913 1913]'; 
tformPivcc = maketform('projective',[u v],[x y]);

u = [397 1636 424 1612]'; v = [1073 1073 1964 1964]'; x = [231 1814 231 1814]';  y =[753 753 1913+63 1913+63]'; 
tformPivsurfcc = maketform('projective',[u v],[x y]);

% imgPivcc_rect=imtransform(pvcc,tformPivcc,'XYScale',1);
imgPivcc_rect=imtransform(imgPivcc,tformPivcc,'XYScale',1);
imgPivcc_rect_cut=imgPivcc_rect(1:2086,62:2015);

% IM1 = imgPivcc_rect_cut - medfilt2(imgPivcc_rect_cut,[5 5]);

%% Pivsurfcc
% imgPivsurfcc_rect=imtransform(pscc,tformPivsurfcc,'XYScale',1);
imgPivsurfcc_rect=imtransform(imgPivsurfcc,tformPivsurfcc,'XYScale',1);
imgPivsurfcc_rect_cut=imgPivsurfcc_rect(542:end-63,352:2305);


figure, imagesc(imgPivsurfcc_rect_cut+imgPivcc_rect_cut), colormap(bone), caxis([0 2500])
figure, imagesc(imgPivsurfcc_rect_cut), colormap(bone), caxis([0 2500])
figure, imagesc(imgPivcc_rect_cut), colormap(bone), caxis([0 1000])



















