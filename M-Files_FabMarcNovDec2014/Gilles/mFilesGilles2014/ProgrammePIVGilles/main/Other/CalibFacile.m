clear all
close all


n = 2048;
m = 2048;
rawFrame = 'E:\data\CamFinale\Calibration\Movie5_Scene1_pivcc_0.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
Pivcc = matFrame.img;

rawFrame = 'E:\data\CamFinale\Calibration\Movie5_Scene1_pivsurfcc_0.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
Pivsurfcc = matFrame.img;

% 
% figure, imagesc(Pivcc), colormap(bone), caxis([0 500])
% figure, imagesc(rot90(Pivsurfcc,0)), colormap(bone), caxis([0 3000])
% 

% u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1088 1091 570]'; 
% u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1088 1091 510]'; 
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]'; 

y=y-10;
% tform = maketform('affine',[u v],[x y]);
tform = maketform('affine',[x y],[u v]);
imgPivsufcc=imtransform(Pivsurfcc,tform,'Xdata',[1 2048],'Ydata',[1 2048], 'FillValues',-1);

% figure, imagesc(imgPivsufcc), colormap(bone), caxis([0 3000])

imgTest(1:763,:)=imgPivsufcc(1:763,:);
imgTest(764:2048,:)=Pivcc(764:2048,:)*3;
figure, imagesc(imgTest), colormap(bone), caxis([0 2000])