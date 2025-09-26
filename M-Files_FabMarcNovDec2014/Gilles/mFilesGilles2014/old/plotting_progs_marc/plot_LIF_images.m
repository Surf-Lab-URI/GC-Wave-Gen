%% Comparateur de resultats
clear all
close all
 for image_pair_number = 360:390
num_of_digits=4

load(['F:\data\LIF\ExpLC3_dt20ms_cc_LIF_1\RawImages\Pivcc\ExpLC3_dt20ms_cc_LIF_1_Pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a'])

figure, imagesc(imgPivcc), colormap(bone), caxis([0 1500])


 end