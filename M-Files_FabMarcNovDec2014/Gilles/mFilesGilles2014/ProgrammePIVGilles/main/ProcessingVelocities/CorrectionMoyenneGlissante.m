clear all
close all
num_of_digits = 4;
path = ['E:\ComputedVelocities\20140804\LC1_1\ExpLC1_1'];
start_image_pair_number = 363;
end_image_pair_number = 371;
GS=4;
x=[GS:GS:2048-GS]; y=x;
mg=4;

for image_pair_number = start_image_pair_number:end_image_pair_number
%     image_pair_number = 365;
    for i=image_pair_number-mg:image_pair_number+mg
        load([path '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'],i) '.mat']);
        delta_z_interp = interpTransfoInverseNew(compVel1.delta_z,SUMASK);
        delta_z_moy(:,i-image_pair_number+mg+1) = nanmean(delta_z_interp,2);
    end
    delta_z_corr = nanmean(delta_z_moy,2);
    delta_z_corr_repmat=nan(size(SUMASK,2));
    for j=1:size(SUMASK,2)
        delta_z_corr_repmat(:,j)=delta_z_corr;
    end
    delta_z_corr_interp = interpTransfoNew(delta_z_corr_repmat,SUMASK);
    load([path '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) '.mat']);
    delta_zcorr = compVel1.delta_z-delta_z_corr_interp;
    
    %% Visualisation
%     figure, imagesc(compVel1.delta_z.*compVel1.mask), colorbar, caxis([-5 5])
%     hold on, quiver(compVel1.delta_x.*compVel1.mask,compVel1.delta_z.*compVel1.mask)
    figure, imagesc(delta_zcorr), colorbar, caxis([-5 5])
    hold on, quiver(compVel1.delta_x.*compVel1.mask,delta_zcorr)
    
end

compVel1.delta_z=delta_zcorr;

