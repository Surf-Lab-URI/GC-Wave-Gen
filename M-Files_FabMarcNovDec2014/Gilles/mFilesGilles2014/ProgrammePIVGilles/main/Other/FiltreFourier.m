tic,
clear all
close all


%% Parameter
exp_name = 'Test';

%%
num_of_digits = 4;

IntrWndw = [128 64 32 16 8]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8 4]; 

save_path = ['E:\ComputedVelocities\Movie7_Scene2\'];
start_image_pair_number = 400;
end_image_pair_number = 450;



%%
% for image_pair_number = start_image_pair_number:end_image_pair_number
    image_pair_number =400;
    n = 2048;
    m = 2048;
    rawFrame = ['E:\Movie7_Scene2\Movie7_Scene2_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    IM1 = matFrame.img;

    for i=1:size(IM1,1)
        F(i,:)=fft(IM1(i,:));
        F(i,2:15)=0;
        IM1f(i,:)=real(ifft(F(i,:)));
    end
    
%     figure, plot(IM1(i,:))
%     figure, plot(IM1f(i,:))
    
    figure, imagesc(IM1), colormap(bone), caxis([0 500])
    figure, imagesc(IM1f), colormap(bone), caxis([0 500])

    
% end