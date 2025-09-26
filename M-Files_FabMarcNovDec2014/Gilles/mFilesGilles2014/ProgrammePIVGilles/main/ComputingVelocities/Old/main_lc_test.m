tic
%%
n = 2048;
m = 2048;
rawFrame = 'E:\Movie7_Scene2\Movie7_Scene2_pivcc_0760.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM1 = matFrame.img;
rawFrame = 'E:\Movie7_Scene2\Movie7_Scene2_pivcc_0761.raw';
matFrame = saveToMatSingleFrame(rawFrame,n,m);
IM2 = matFrame.img;
%%
IntrWndw = [128 64 32 16 8]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8 4]; 
%%
    mask1 = nan(size(IM1));
    for i=900:size(IM1,1)
        mask1(i,:) = 1;
    end
compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask1, IntrWndw, GrdSpc);
