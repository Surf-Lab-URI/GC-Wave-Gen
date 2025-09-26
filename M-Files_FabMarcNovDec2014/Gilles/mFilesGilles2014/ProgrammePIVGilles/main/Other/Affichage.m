clear all
close all
num_of_digits = 4;
path = ['E:\data\20140826\LC1_4\Movie10_Scene7'];
start_image_pair_number = 363;
end_image_pair_number = 366;
GS=4;
x=[GS:GS:2048-GS]; y=x;

for image_pair_number = start_image_pair_number:end_image_pair_number
    load(['E:\ComputedVelocities\20140826\LC1_4\ExpLC1_4_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat']);
    
    
        %% Comparaison compVel vs compVel1
        image_letter='a';
        rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2) '.raw'];
        matFrame = saveToMatSingleFrame(rawFrame,2048,2048);
        IM1 = matFrame.img;
    
        image_letter='b';
        rawFrame = [path '_pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number*2+1) '.raw'];
        matFrame = saveToMatSingleFrame(rawFrame,2048,2048);
        IM2 = matFrame.img;
    
        figure, imagesc(IM1), colormap(bone), caxis([0 500])
%         hold on, quiver(x,y,compVel.delx_ints.*compVel.MASK,compVel.dely_ints.*compVel.MASK)
        hold on, quiver(x,y,compVel1.delta_x.*compVel1.mask,compVel1.delta_z.*compVel1.mask,'r')
%         hold on, quiver(x,y,delx.*mask3(4:4:end-4,4:4:end-4),dely.*mask3(4:4:end-4,4:4:end-4),'r')
        figure, imagesc(IM2), colormap(bone), caxis([0 500])
    
    
%     %% Recherche LCs sans correction
%     figure, imagesc(compVel.dely_ints.*compVel.MASK), colorbar, caxis([-5 5])
%     hold on, quiver(compVel.delx_ints.*compVel.MASK,compVel.dely_ints.*compVel.MASK)
%     
%     figure, imagesc(compVel1.delta_z.*compVel1.mask), colorbar, caxis([-5 5])
%     hold on, quiver(compVel1.delta_x.*compVel1.mask,compVel1.delta_z.*compVel1.mask,'r')
%     
%     %% Recherche LCs avec correction
%     figure, imagesc((compVel.dely_ints+dSUInterp).*compVel.MASK), colorbar, caxis([-5 5])
%     hold on, quiver(compVel.delx_ints.*compVel.MASK,(compVel.dely_ints+dSUInterp).*compVel.MASK)
%     
%     figure, imagesc((compVel1.delta_z+dSUInterp).*compVel1.mask), colorbar, caxis([-5 5])
%     hold on, quiver(compVel1.delta_x.*compVel1.mask,(compVel1.delta_z+dSUInterp).*compVel1.mask,'r')
end