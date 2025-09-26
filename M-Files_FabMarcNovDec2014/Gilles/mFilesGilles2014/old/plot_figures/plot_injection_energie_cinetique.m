%% Export figures - Injection d'energie cinetique


tic,
clear all
close all


%% Parameter
exp_name = 'LC1_dt25ms_1';
deltaT = 25d-3; %s
path = '\\beo\data\';
save_path = 'F:\figures\';
num_of_digits = 4;
%% Cameras resolution
piv_res = 40d-6; %m/pix
vec_res = 4 * piv_res;
pivsurf_res = 102d-6; %m/pix
piv_delta_t = 1/7.2; %sec
pivResReal = 40d-6; %m/pix
pivRes = 1;%40d-6; %m/pixel
iws = 16; % initial widow size
vResLfv =  1.3543e-04/pivResReal; %pixel de piv/pixel de lfv;  %lfv resolution
hResLfv = 1.2912e-04/pivResReal; %pixel de piv/pixel de lfv

psSWL =  (2048 - 575 + 1) * pivResReal;
XX = 0:pivResReal:3785*pivResReal-pivResReal;
YYL = 2048*pivResReal-pivResReal;
% YY = - psSWL : pivResReal : YYL - psSWL;
YY = -YYL:pivResReal:0;

%%

ImgDebut=300;
ImgFin=600;

U2=zeros(255,ImgFin-ImgDebut+1);
W2=zeros(255,ImgFin-ImgDebut+1);
XXtemps=0;      


for image_number = ImgDebut:ImgFin
    load([path 'Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number)]);
    %% Interpolation de l'image
    delx=0;
    delxOr=0;
    surface=0;
    compVel.delx_ints=compVel.delx_ints*pivResReal;
    compVel.dely_ints=compVel.dely_ints*pivResReal;
    compVel.delxOrb=compVel.delxOrb*pivResReal;
    compVel.delyOrb=compVel.delyOrb*pivResReal;
    XXtemps(image_number-ImgDebut+1)=image_number*(1/7.2);
    for col=1:size(compVel.delx_int,2)
        xi=0;
        xt=0;
        surface(col)=find(compVel.MASK(:,col)==1,1);
        xt=surface(col):1:size(compVel.delx_int,1);
        delx(1:length(xt),col)=compVel.delx_ints(surface(col):end,col);
        dely(1:length(xt),col)=compVel.dely_ints(surface(col):end,col);
        delxOr(1:length(xt),col)=compVel.delxOrb(surface(col):end,col);
        delyOr(1:length(xt),col)=compVel.delyOrb(surface(col):end,col);
    end
    %     U(1:size(delx,1),image_number-ImgDebut+1)=nanmean(delx,2);
    %     UOrb(1:size(delx,1),image_number-ImgDebut+1)=nanmean(delx-delxOr,2);
    U2(1:size(delx,1),image_number-ImgDebut+1)=nanmean((delx-delxOr).*(delx-delxOr),2);
    W2(1:size(dely,1),image_number-ImgDebut+1)=nanmean((dely-delyOr).*(dely-delyOr),2);
end
E=U2+W2;

%% Axes dans la bonne unite
XXvel = XX(8:8:end-8);
YYvel = YY(8:8:end-8);


%% horizontal
figure, imagesc(XXtemps,YYvel,flipud(E))
caxis([0 1e-6])


col = jet(64); tmp = linspace(0,1,64)';
for n = 1:3, col(:,n) = interp1( 10.^tmp, col(:,n), 1+9*tmp, 'linear'); end
colormap(col)


% colorbar
h = colorbar;

% cmap_u = [floor(10*nanmin(u(:)))/10 ceil(10*nanmax(u(:)))/10];
% set(h,'YTick',[0 0.5 1]*1e-6);
% set(get(h,'ylabel'),'Interpreter','latex','String', '{$$u~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% set(h,'ylabel','E','FontSize',14, 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('t(s)', 'FontSize',14);
ylabel('z(m)', 'FontSize',14)
% background color and size
set(gcf, 'color', 'w');
set(gcf, 'OuterPosition', [399-150   242-150   568+240   502])

% strr = ['export_fig ' save_path 'Exp' exp_name '_Ec_'' -jpg -nocrop -r500 -q100'];
eval(['export_fig ' save_path 'Exp' exp_name '_Ec_moinsOrb -jpg -nocrop -r500 -q100'])

% eval(['export_fig ' save_path 'Exp' exp_name '_difference_w_Orbw_' ' -jpg -nocrop -r500 -q100'])
toc