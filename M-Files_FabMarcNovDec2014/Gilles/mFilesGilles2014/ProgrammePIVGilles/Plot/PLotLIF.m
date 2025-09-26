
path = ['E:\data\CamFinale\PivLif\PivLif1\Movie7_Scene1']
path = ['E:\data\CamFinale\PivLif\PivLif2\Movie7_Scene2']
path = ['E:\data\CamFinale\PivLif\PivLif3\Movie7_Scene3']
path = ['E:\data\CamFinale\PivLif\PivLif4\Movie7_Scene7']
num_of_digits = 4;
n=2048; m=2048;
ResPivReal = 6.4101e-05;
inter=[400:1200];
xx=[1:2048]*ResPivReal;
yy=(inter-720)*ResPivReal;


image_pair_number=353;
image_pair_number/7.2-21
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')

image_pair_number=356;
image_pair_number/7.2-21
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')

image_pair_number=359;
image_pair_number/7.2-21
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')

image_pair_number=362;
image_pair_number/7.2-21
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')

image_pair_number=365;
image_pair_number/7.2-21
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')

image_pair_number=368;
rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
figure, imagesc(xx,yy,IM1(inter,:)), colormap(bone), caxis([0 500]), axis('equal')
