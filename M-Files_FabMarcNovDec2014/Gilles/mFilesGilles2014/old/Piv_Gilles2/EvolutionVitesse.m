%% Calcul de l'evolution du profil des vitesses en fonction du temps
clear all
% close all

ImgDebut = 300;
ImgFin = 600;

exp_name = 'LC2_dt25ms_1';
deltaT = 25d-3; %s

num_of_digits = 4;
path = ['\\beo\data\Exp' exp_name '\ComputedVelocities\'];
save_path = 'F:\figures\';



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







U=zeros(255,ImgFin-ImgDebut+1);

for image_number = ImgDebut:ImgFin
    load([path 'Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_number)]);
    %% Interpolation de l'image
    delx=0;
    delxOr=0;
    surface=0;
     compVel.delx_ints=compVel.delx_ints*pivResReal/deltaT;
    
    for col=1:size(compVel.delx_int,2)
        xi=0;
        xt=0;
        surface(col)=find(compVel.MASK(:,col)==1,1);
        xt=surface(col):1:size(compVel.delx_int,1);
        delx(1:length(xt),col)=compVel.delx_ints(surface(col):end,col);
%         dely(1:length(xt),col)=compVel.dely_ints(surface(col):end,col);
        delxOr(1:length(xt),col)=compVel.delxOrb(surface(col):end,col);
%         delyOr(1:length(xt),col)=compVel.delyOrb(surface(col):end,col);
    end
    U(1:size(delx,1),image_number-ImgDebut+1)=nanmean(delx,2);
    UOrb(1:size(delx,1),image_number-ImgDebut+1)=nanmean(delx-delxOr,2);
%     U2(1:size(delx,1),image_number-ImgDebut+1)=nanmean(delx.*delx,2);
%     W2(1:size(dely,1),image_number-ImgDebut+1)=nanmean(dely.*dely,2);
end
% E=U2+W2;

%%
% figure, imagesc(U)
% for image_number = ImgDebut:ImgFin
%     plot(U(:,image_number-ImgDebut+1));
%     M(image_number)=getframe;
% end
% movie(M,1,1)
% % 
% 
% % %%

XXvel = XX(8:8:end-8);
YYvel = YY(8:8:end-8);

%LC3
figure, plot(U(:,2))
hold on, plot(U(:,70),'g'),plot(U(:,80),'c'), plot(U(:,90),'r'), plot(U(:,100),'k'), plot(U(:,103),'m')




% cmap_u = [floor(10*nanmin(u(:)))/10 ceil(10*nanmax(u(:)))/10];
% set(h,'YTick',[0 0.5 1]*1e-6);
% set(get(h,'ylabel'),'Interpreter','latex','String', '{$$u~(m/s)$$}','FontSize',14, 'Position',get(get(h,'ylabel'),'Position') - [0 .01 0] , 'rotation',90);
% set(h,'ylabel','E','FontSize',14, 'rotation',90);
% x and z labels
% hold on, plot(XX,(size(surf1.img,1) - surf1.z_s(end:-1:1)+1)*pivResReal - psSWL,'r','LineWidth',2)
set(gca,'xdir','normal'); set(gca,'ydir','normal');

% x and z labels
xlabel('z(m)', 'FontSize',14);
ylabel('U(m/s)', 'FontSize',14);
% background color and size
set(gcf, 'color', 'w');


% strr = ['export_fig ' save_path 'Exp' exp_name '_Ec_'' -jpg -nocrop -r500 -q100'];
eval(['export_fig ' save_path 'Exp' exp_name 'lfv -jpg -nocrop -r500 -q100'])
