%% Synthetic 
clear
clc
close all

load('SyntheticSurf_Crapper_Noisier1.0cm')

figure(1)
imagesc(im)
hold on
plot(x/mpp+1,v_offset - y/mpp,'-r');
[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface] = FindWaterSurface(im);
figure(1)
hold on
plot(XPIVSurfW1_Surface, PIVSurfW1_Surface,'-r', 'LineWidth',2);

set(gca,'DataAspectRatio',[1 1 1])
% ylim([1000,1500])

%%
clear
clc
close all

load('CapTestImg1094.mat')
figure(1)
imagesc(PIVSurfW1_CamAngle,[0,70])
hold on

[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface] = FindWaterSurface(PIVSurfW1_CamAngle);

figure(1)
hold on
plot(XPIVSurfW1_Surface, PIVSurfW1_Surface,'-r', 'LineWidth',2);

%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle)
    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = CrapperOptimized_FindSurface(PIVSurfW_CamAngle(1:2800,X), 50, 1);
    PIVSurf_Surface_Raw = imSurf.surface;
    f = gcf;
    figure(1)
    hold on
    plot(X,PIVSurf_Surface_Raw,'-g')
    
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
    plot(X,PIVSurf_Surface_Int,'-k')
    % PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
    % plot(X,PIVSurf_Surface_Int,'-m')
    figure(2)
    plot(500*diff(PIVSurf_Surface_Int,2))
    [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 2d3);
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