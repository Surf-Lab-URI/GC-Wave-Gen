tic,
clear all
close all


%% Parameter
exp_name = 'Test';

%%
num_of_digits = 4;

IntrWndw = [128 64 32 16 8]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8 4]; 

save_path = ['E:\ComputedVelocities\Movie4_Scene2\'];
start_image_pair_number = 390;
end_image_pair_number = 409;



%%
for image_pair_number = start_image_pair_number:end_image_pair_number
    
    n = 2048;
    m = 2048;
    rawFrame = ['E:\Movie4_Scene2\Movie4_Scene2_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM1 = matFrame.img;
    rawFrame = ['E:\Movie4_Scene2\Movie4_Scene2_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM2 = matFrame.img;
    
    mask1 = nan(size(IM1));
    for i=730:size(IM1,1)
        mask1(i,:) = 1;
    end
    
    compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask1, IntrWndw, GrdSpc);
    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
 end
toc