clear
clc
close all

mpp = 41.782d-6/0.67; % meters per pixel

rows = 3113; 
cols = 4176;

imsurf = zeros(rows,cols);

lambda = 1e-2;
a = lambda*0.73;

A = 2*lambda/pi/a*((1+pi^2*a^2/4/lambda^2)^0.5-1);

alpha = 0:0.0002:ceil(cols*mpp/lambda);

x = lambda*(alpha-2/pi*A*sin(2*pi*alpha)./(1+A^2+2*A*cos(2*pi*alpha)));
y = lambda*(2/pi-2/pi*A*cos(2*pi*alpha)./(1+A^2+2*A*cos(2*pi*alpha)));

plot(x/mpp,y/mpp)
axis equal

v_offset = 1500;
figure
imax = find(round(x/mpp)+1==cols);
is = sub2ind([rows,cols],v_offset-round(y(1:imax(1))/mpp),round(x(1:imax(1))/mpp)+1);
imsurf(is)=1;

imbinfill = imfill(logical(imsurf),sub2ind([rows,cols],rows,cols));
imagesc(imbinfill)

figure
im = zeros(rows,cols);
im(imbinfill) = 60;
im(~imbinfill) = 5;
noise = randn(rows,cols);
im = max(im+20*noise,0);
im = imgaussfilt(im,5);
imagesc(im)
s = sprintf("SyntheticSurf_Crapper_Noisier%.1fcm.mat", 100*lambda);
save(s)
%%
[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface] = FindWaterSurface(im);
hold on
plot(XPIVSurfW1_Surface, PIVSurfW1_Surface,'-r');

set(gca,'DataAspectRatio',[1 1 1])
ylim([1000,1500])

%%
function mat = gauss2d(mat, sigma, center)
    gsize = size(mat);
    for r=1:gsize(1)
        for c=1:gsize(2)
            mat(r,c) = gaussC(r,c, sigma, center);
        end
    end
end

function val = gaussC(x, y, sigma, center)
    xc = center(1);
    yc = center(2);
    exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
    val       = (exp(-exponent));
end

function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle)
    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = CrapperOptimized_FindSurface(PIVSurfW_CamAngle(1:2800,X), 3, 5);
    PIVSurf_Surface_Raw = imSurf.surface;
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
    [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1d5);
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