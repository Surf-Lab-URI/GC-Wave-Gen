close all
rawFrame_pivcc = '\\beo\data\Calib_cc\Movie13_Scene22_pivcc_0.raw';
rawFrame_pivsurfcc = '\\beo\data\Calib_cc\Movie13_Scene22_pivsurfcc_0.raw';
n = 2048;
m = 2048;
p = saveToMatSingleFrame(rawFrame_pivcc,n,m);
ps = saveToMatSingleFrame(rawFrame_pivsurfcc,n,m);

pivcc = p.img;
pivsurfcc = rot90(ps.img,2);

figure, imagesc(pivcc), colormap(bone), caxis([0 500])
figure, imagesc(pivsurfcc), colormap(bone), caxis([0 1500])

