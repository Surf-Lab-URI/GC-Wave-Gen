% Function findSurfaceLfvFabMarc
%% Object:
% detects surface on LFV LIF images
%% Arguments: 
% imgLfv (2D double)
%% Result: 
% imSurf.x_s (x-axis vector in pixels), imSurf.z_s (surface elevation vector in pixels
% (upside down, upright surface is 2048 - imSurf.z_s + 1))
%% Author:
% Fabrice Veron, modified by Marc Buckley
%% Last update:
% 06/26/2013
%% Example:
% load('E:\data\Exp4\RawImages\Lfv\Exp4_Lfv_1004.mat')
% imSurf = findSurfaceLfvFabMarc(imgLfv);
% figure, imagesc(imgLfv), colormap(bone), hold on, plot(imSurf.x_s, imSurf.z_s, 'r')
%%
function imSurf = findSurfaceLfvFabMarc(imgLfv,thresh)

badIm = 0;
img = imgLfv;
% [~,thresh] = edge(img,'prewitt');
bw = edge(img,'prewitt',thresh);
% bw = edge(img,'prewitt');
bw = double(bw);
surface = nan(size(img,2));
for i=1:size(img,2)
    surfTemp = find(bw(:,i),1,'first');
    if isempty(surfTemp)
        surface(i) = nan;
        badIm = badIm+1;
    else
        surface (i) = surfTemp;
    end
end

s=(despike_fab(surface));
smth_surface= smoothn(s,100);
% h=fdesign.lowpass('Fp,Fst,Ap,Ast',0.08,0.5,1,100);
%     d=design(h,'equiripple'); %Lowpass FIR filter
%     y=filtfilt(d.Numerator,1,s); %zero-phase filtering
% smth_surface = s;

imSurf.x_s = 1:size(img,2);
imSurf.z_s = smth_surface;
imSurf.z_s_raw = surface;
imSurf.badImColCounter = badIm;