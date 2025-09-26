tic,
clear all
close all


%% Parameter
exp_name = 'LC3_dt25ms_1';
deltaT = 25d-3; %s
path = '\\beo\data\';
save_path = 'F:\figures\';
num_of_digits = 4;
piv_res = 40d-6; %m/pix
vec_res = 4 * piv_res;
pivsurf_res = 102d-6; %m/pix
piv_delta_t = 1/7.2; %sec

pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
iws = 16; % initial widow size
vResLfv =  1.3543e-04/pivResReal; %pixel de piv/pixel de lfv;  %lfv resolution
hResLfv = 1.2912e-04/pivResReal; %pixel de piv/pixel de lfv

image_pair_number = 790


psSWL =  (2048 - 575 + 1) * pivResReal;
XX = 0:pivResReal:3785*pivResReal-pivResReal;
YYL = 2048*pivResReal-pivResReal;
YY = - psSWL : pivResReal : YYL - psSWL;

% YY = fliplr(YY);

%% plot velocity field

load(['F:\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat']);

XXvel = XX(8:8:end-8);
YYvel = YY(8:8:end-8);



%% horizontal
u = compVel.delx_int.*compVel.MASK .* piv_res / deltaT;
figure, imagesc(XXvel,YYvel,flipud(u))
caxis([-0.17 0.17])
% colorbar
h = colorbar;
% cmap_u = [floor(10*nanmin(u(:)))/10 ceil(10*nanmax(u(:)))/10];
% set(h,'YTick',[cmap_u(1) 0 1 2 cmap_u(2)]);
set(get(h,'ylabel'),'Interpreter','latex','String', '{$$u~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

eval(['export_fig ' save_path 'Exp' exp_name '_u_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])

%% vertical
w = -compVel.dely_int.*compVel.MASK .* piv_res / deltaT;
figure, imagesc(XXvel,YYvel,flipud(w))
caxis([-0.17 0.17])
% colorbar
h = colorbar;
% cmap_w = [floor(10*nanmin(w(:)))/10 ceil(10*nanmax(w(:)))/10];
% set(h,'YTick',[cmap_u(1) 0 1 2 cmap_u(2)]);
set(get(h,'ylabel'),'Interpreter','latex','String', '{$$w~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

eval(['export_fig ' save_path 'Exp' exp_name '_w_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])

%% horizontal vitesse-orbitales

u = (compVel.delx_int-compVel.delxOrb).*compVel.MASK .* piv_res / deltaT;
figure, imagesc(XXvel,YYvel,flipud(u))
caxis([-0.1 0.1])
% colorbar
h = colorbar;
% cmap_u = [floor(10*nanmin(u(:)))/10 ceil(10*nanmax(u(:)))/10];
% set(h,'YTick',[cmap_u(1) 0 1 2 cmap_u(2)]);
set(get(h,'ylabel'),'Interpreter','latex','String', '{$$u^{\prime}~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

eval(['export_fig ' save_path 'Exp' exp_name '_difference_u_Orbu_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])
%%
w = -(compVel.dely_int-compVel.delyOrb).*compVel.MASK .* piv_res / deltaT;
% figure, imagesc(XXvel,YYvel,flipud(w),'AlphaData',~isnan(flipud(w)),'AlphaDataMapping','none') 
figure, imagesc(XXvel,YYvel,flipud(w)),
caxis([-0.1 0.1])
% colorbar
h = colorbar;
% cmap_w = [floor(10*nanmin(w(:)))/10 ceil(10*nanmax(w(:)))/10];
% set(h,'YTick',[cmap_u(1) 0 1 2 cmap_u(2)]);
set(get(h,'ylabel'),'Interpreter','latex','String', '{$$w^{\prime}~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('X(m)', 'FontSize',14);
ylabel('Z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])
eval(['export_fig ' save_path 'Exp' exp_name '_difference_w_Orbw_' sprintf(['%0' num2str(num_of_digits) 'd'],image_pair_number) ' -jpg -nocrop -r500 -q100'])


% map = colormap('jet');
% map1 = [map; [1 1 1]];
% ww = w;
% ww(isnan(w)) = 1000;


