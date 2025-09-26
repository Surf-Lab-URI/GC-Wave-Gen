%% Comparateur de resultats
clear all
close all
exp_name = 'LC3';
 for image_pair_number = 370
num_of_digits = 4

load(['F:\data\ExpLC3_dt20ms_cc_1\ComputedVelocities\ExpLC3_dt20ms_cc_1_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt'])

% load(['F:\data\ExpLC3_dt20ms_cc_2\ComputedVelocities\ExpLC3_dt20ms_cc_2_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt'])
% langmuir on 377 maybe

% load(['F:\data\ExpLC3_dt20ms_cc_1\ComputedVelocities\ExpLC3_dt20ms_cc_1_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt'])
% langmuir on 371 maybe
% set(gca,'xdir','normal'); set(gca,'ydir','normal');
pivcc_res = 130e-6 * 100; % cm/pixel
vec_res = pivcc_res * 4;
% figure, imagesc(compVel.dely_ints.*compVel.MASK),colorbar
% figure, imagesc(curl(compVel.delx_ints.*compVel.MASK,compVel.dely_ints.*compVel.MASK)),colorbar, caxis([-max(abs(caxis)) max(abs(caxis))])
% figure, imagesc(compVel.dcor.*compVel.MASK)
xx = compVel.delx_ints.*compVel.MASK.*200e-6/(20e-3)*5;  %results in 5*m/s
yy = compVel.dely_ints.*compVel.MASK.*200e-6/(20e-3)*5;
xx_sub = xx(180+17:180+44,100+28:100+68);
% yy_sub = yy(180:240,100:200);
yy_sub = yy(180+17:180+44,100+28:100+68);
% figure, imagesc(xx), colorbar
[X,Y] = meshgrid(0:vec_res:vec_res*(41-1),-1.14:vec_res:vec_res*(6));
figure, quiver(X,Y,flipud(xx_sub), flipud(-yy_sub),0)

hold on, quiver(0.5,0.2,0.01*5,0,0)  % cm/s

% bb = [27.5 67.6 19.7 43.9];
axis([0 2.1 -1.1 0.3]);
% axis(bb)

% set(gca,'xdir','normal'); set(gca,'ydir','normal');
% pivcc_res = 130e-6 * 100; % cm/pixel
% vec_res = pivcc_res * 4;

% set(gca,'YTick',[0 5*vec_res 2*5*vec_res 3*5*vec_res 4*5*vec_res]);
% x and z labels
xlabel('X(cm)', 'FontSize',14);
ylabel('Z(cm)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502+240])
save_path = 'F:\mFiles\ProgrammePIVGilles\';
% eval(['export_fig ' save_path 'Exp' exp_name '_quiver_LC_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])



% 
% dcor = compVel.dcor;
% MASK = compVel.MASK;
% dcor1 = 50 * ones(size(dcor)); % matrix of 50s
% dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
% % figure, imagesc(dcor1), caxis([0 1]), colorbar
% numel(find(dcor<0.5))/numel(find(dcor1<40))


% load(['F:\data\ExpLC2_dt25ms_cc_1\ComputedVelocities\ExpLC2_dt25ms_cc_1_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)])
% 
% figure, imagesc(compVel.delx_ints.*compVel.MASK),colorbar
% figure, imagesc(compVel.dely_ints.*compVel.MASK ),colorbar%, caxis([-3 3])
% figure, imagesc(compVel.dcor.*compVel.MASK)
% figure, quiver(compVel.delx_ints.*compVel.MASK,(compVel.dely_ints- mean(compVel.dely_ints(:))).*compVel.MASK,2)
% 
% 
% dcor = compVel.dcor;
% MASK = compVel.MASK;
% dcor1 = 50 * ones(size(dcor)); % matrix of 50s
% dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
% % figure, imagesc(dcor1), caxis([0 1]), colorbar
% numel(find(dcor<0.5))/numel(find(dcor1<40))

 end