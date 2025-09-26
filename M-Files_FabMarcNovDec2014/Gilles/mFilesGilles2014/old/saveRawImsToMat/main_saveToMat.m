n = 2048;
m = 2048;
rawFrame = 'F:\silver_part_piv\M1S4\Movie1_Scene4_piv1_0700.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;

rawFrame = 'F:\silver_part_piv\M1S4\Movie1_Scene4_piv1_0701.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM2 = matFrame.img;