%% Synthetic 
clear
clc
% close all

load('SyntheticSurf_Crapper_MultiscaleNoise1.0cm.mat')



figure(1)
hold off
imagesc(im,[0,100])
hold on
p1 = plot(x/mpp+1,v_offset - y/mpp,'-b','DisplayName','Crapper Capillary Wave','LineWidth',1);
[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface,p3] = FindWaterSurface(im);
figure(1)
hold on
p2 = plot(XPIVSurfW1_Surface, PIVSurfW1_Surface,'-r', 'LineWidth',2,'DisplayName','Surface Detection Output');
set(gca,'DataAspectRatio',[1 1 1])
axis off
set(gca,'FontSize',24)
% title(s,'Interpreter','none')

xl = xlim;
yl = ylim;

lsbm = 1e-2; % length of scale bar in meters
lsb = lsbm/mpp;
xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
ysb = (yl(2) - (yl(2)-yl(1))*0.05)*[1 1];
try
    delete(sb)
end
sb = plot(xsb,ysb,'-k', 'LineWidth',10);


sbl = sprintf('%d cm',lsbm*100);
text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex')


legend([p1, p2, p3], 'Interpreter','latex')


% figure(3)
% plot(-diff(PIVSurfW1_Surface,2)*100+mean(PIVSurfW1_Surface),'-m')
% ylim([1000,1500])

%% Real
clear
clc
% close all
% mpp = 6.236119402985075e-5;
load('CapTestImg1094.mat')


figure(2)
hold off
imagesc(PIVSurfW1_CamAngle,[0,100])
hold on

[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, p2] = FindWaterSurface(PIVSurfW1_CamAngle);

figure(2)
hold on
p1 = plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, '-r', 'LineWidth',2,'DisplayName','Surface Detection Output');
set(gca,'DataAspectRatio',[1 1 1])
axis off
set(gca,'FontSize',24)
% title(s,'Interpreter','none')

ylim([1750,2150])

xl = xlim;
yl = ylim;

lsbm = 1e-2; % length of scale bar in meters
lsb = lsbm/mpp;
xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
ysb = (yl(2) - (yl(2)-yl(1))*0.15)*[1 1];
try
    delete(sb)
    delete(sbt)
end
sb = plot(xsb,ysb,'-k', 'LineWidth',10);


sbl = sprintf('%d cm',lsbm*100);
sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex');

mplot = plot([XPIVSurfW1_Surface(1),XPIVSurfW1_Surface(end)],[mean(PIVSurfW1_Surface),mean(PIVSurfW1_Surface)],'--m','LineWidth',2, 'DisplayName','Mean Water Level');


delete(p2)


legend([p1, mplot], 'Interpreter','latex')

% figure(3)
% plot(-diff(PIVSurfW1_Surface,2)*10000+2000,'-m')
%% Slope and Length Scales
xl = round(xlim);

ixi = find(XPIVSurfW1_Surface == xl(1));
ixf = find(XPIVSurfW1_Surface == xl(2));

ysurf = PIVSurfW1_Surface(ixi:ixf);
xsurf = XPIVSurfW1_Surface(ixi:ixf);

eta_pix = ysurf-mean(PIVSurfW1_Surface);

[minlocs,p] = islocalmax(ysurf);

[m,ix0] = max(p);

x0 = xsurf(ix0);
plot(x0,mean(ysurf),'.r','MarkerSize',15)

L = sqrt(sum(eta_pix.^2.*(xsurf-x0).^2)/sum(eta_pix.^2))

% plot([x0-L/2, x0+L/2],[mean(ysurf-20),mean(ysurf-20)],'LineWidth',6,'DisplayName','$L$')
plot([xl(1),xl(2)],[mean(PIVSurfW1_Surface),mean(PIVSurfW1_Surface)],'--m','LineWidth',4)

eta_x = diff(eta_pix);
[eta_x_max,iexm] = max(eta_x);
plot(xsurf(iexm),ysurf(iexm),'.r','MarkerSize',15)
plot(xsurf(iexm)+[-50,50],ysurf(iexm)+[-50,50]*eta_x_max,'--b','LineWidth',4)
plot(xsurf(iexm)+[15,-15,-15],ysurf(iexm)+[15,15,-15]*eta_x_max,'-b','LineWidth',4)
s = sprintf('%.2f',eta_x_max);
text(xsurf(iexm)-15-3,ysurf(iexm),s,'FontSize',20,'HorizontalAlignment','right','Interpreter','latex')
ylim([1950,2100])
legend([p1, p2, mplot], 'Interpreter','latex')

%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface,p3] = FindWaterSurface(PIVSurfW_CamAngle)
    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = CrapperOptimized_FindSurface(PIVSurfW_CamAngle(1:2800,X), [50,40,30,20,10,8,6],[50,40,30,20,10,8],1);
    PIVSurf_Surface_Raw = imSurf.surface;
    f = gcf;
    figure(f)
    hold on
    p3 = plot(X,PIVSurf_Surface_Raw,'-k','LineWidth',5, 'DisplayName','Local Gradient Maximum');
    
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
    % plot(X,PIVSurf_Surface_Int,'-k')

    % PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
    % plot(X,PIVSurf_Surface_Int,'-m')
    % figure
    % plot(500*diff(PIVSurf_Surface_Int,2))
    [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 3d2);
    if length(SP.coefs)>2
        PIVSurf_Surface_W = SP.coefs(2:end-1);
    else
        CC = polyfit([X(1),X(end)],[SP.coefs(1) SP.coefs(2)],1);
        PIVSurf_Surface_W = polyval(CC,[X(1):X(end)]);
    end
    Usurf = X;
    Vsurf = PIVSurf_Surface_W;
    
    %%% Check if bad frame
    BadFramePIVSurfW = 0;
    if imSurf.badFrameBool == 1
        BadFramePIVSurfW = 1;
    end
    
    XPIVSurfW_Surface = Usurf;
    PIVSurfW_Surface = Vsurf;
end