clear all
close all
exp_name = 'Test';
start_image_number = 380;
end_image_number = 390;
path = 'E:\ComputedVelocities\Movie4_Scene2\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];

for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    %     load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    load([path 'Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt']);
    %
%     figure, imagesc(compVel.dely_ints.*compVel.MASK), colorbar, caxis([-5 5])
        figure, imagesc(curl(compVel.delx_ints.*compVel.MASK,compVel.dely_ints.*compVel.MASK)),colorbar, caxis([-0.1 0.1])
    hold on, quiver(compVel.delx_ints.*compVel.MASK,compVel.dely_ints.*compVel.MASK)
    
    % Correction 'barbare' de la vitesse verticale induite par la reflexion
    % de la surface (a refaire en moyennant sur plusieurs images successives ou avec une moyenne d'ensemble)
%     for i = 1:size(compVel.dely_ints,1)
%         meany(i)=mean(compVel.dely_ints(i,:).*compVel.MASK(i,:));
%         dely_cor(i,:)=compVel.dely_ints(i,:).*compVel.MASK(i,:)-meany(i);
%     end
%     figure, imagesc(dely_cor), colorbar, caxis([-3 3])
%     figure, quiver(compVel.delx_ints.*compVel.MASK,dely_cor.*compVel.MASK)
end

%
% for i = 1:size(compVel.dely_ints,1)
%     meany(i)=mean(compVel.dely_ints(i,:).*compVel.MASK(i,:));
% end
% figure, plot(meany)
% % 
% % figure, imagesc(compVel.dely_ints.*compVel.MASK), colorbar, caxis([-5 5])