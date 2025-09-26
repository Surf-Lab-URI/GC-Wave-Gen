% plot_velocityFieldSnapshotColormap
%% Object:
% plots instantaneous velocity field from PIV computation, with the water
% image from the corresponding PivSurf image
%% Functions called: 
%% Result: 
% jpg or png
%% Author:
% Marc Buckley
%% Last update:
% 07/02/2013
%%
clear all
close all
load('E:\MFiles\plotFigs\fabmap.mat');
%% user input
exp_name = '4';
image_pair_number = 1004;
%% constants
path = ['E:\data\Exp' exp_name '\'];
save_path = ['E:\figuresPaper1\'];
num_of_digits = 4;
piv_res = 47.4; %um/pix
vec_res = 4 * piv_res;
pivsurf_res = 99.6; %um/pix
piv_delta_t = 1/7.2; %sec
%% generate Pivsurf image
disp(['pair ' num2str(image_pair_number) ' finding surface...']);
% fuse first frames of pair
image_letter = 'a';
fusedIm = fuseImages(exp_name, [path 'RawImages\Piv1\'],[path 'RawImages\Piv2\'], image_pair_number, num_of_digits, image_letter);
IM1 = fusedIm.fused_im;
% find surface on first fused frame and generate mask for PIV computation
imSurf = findSurface(exp_name, [path 'RawImages\Pivsurf\'], image_pair_number, num_of_digits, image_letter, IM1, piv_res, pivsurf_res);
mask = imSurf.mask;
x_s = imSurf.x_s;
z_s = imSurf.z_s;
disp(['pair ' num2str(image_pair_number) ' surface found']);
%% Processed velocity
% Exp6_compVel_300
procVel = [path 'ProcessedVelocities\Exp' exp_name '_procVel_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) '.mat'];
%
load(procVel)
%
u = PV.U;
mask_inv = nan(size(mask));
mask_inv(isnan(mask)) = 1;
pvs = flipud(imSurf.pivSurf_resized.*mask_inv);
pvs = pvs(PV.z,PV.x);

%% user input
cmp_leng = 64;
cmp_str = ['jet(' num2str(cmp_leng) ')'];
cmp = eval(cmp_str);
cmap_u = [floor(10*nanmin(u(:)))/10 ceil(10*nanmax(u(:)))/10]
% cmap_u = [nanmin(u(:)) nanmax(u(:))]
% cmap_u = [-1 4]
cmap_pvs = [nanmin(pvs(:)) nanmax(pvs(:))]
% cmap_pvs = [0 2000]
%%
% generate rgb masks
mask_air_rgb = repmat(flipud(isnan(mask(PV.z,PV.x))),[1 1 3]);
mask_water_rgb = repmat(flipud(isnan(mask_inv(PV.z,PV.x))),[1 1 3]);
% convert u to rgb
u_gray = mat2gray(u,cmap_u);
u_ind = gray2ind(u_gray,cmp_leng);
u_rgb = ind2rgb(u_ind,cmp) .* mask_water_rgb;
% convert pvs to rgb
pvs_gray = mat2gray(pvs,cmap_pvs);
pvs_ind = gray2ind(pvs_gray,cmp_leng);
pvs_rgb = ind2rgb(pvs_ind,gray(cmp_leng)) .* mask_air_rgb;
% final image
final_im = u_rgb + pvs_rgb;
%% plot and save
figure,imagesc(PV.X,PV.Z,final_im), caxis(cmap_u), colormap(cmp), colorbar
set(gca,'xdir','normal'); set(gca,'ydir','normal');
% colorbar
h = colorbar;
set(h,'YTick',[cmap_u(1) 0 1 2 cmap_u(2)]);
set(get(h,'ylabel'),'Interpreter','latex','String', '{$$u~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% x and z labels
xlabel('X(m)', 'FontSize',14); 
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])
%% save to jpg
vname=@(x) inputname(1);
toto=cmp;
s=vname(cmp)
eval(['export_fig ' save_path 'Exp' exp_name '_PVdotU_' cmp_str '_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])



