% Function findSurface
%% Object:
% detects surface
%% Arguments: 
% img(2D double)
%% Result: 
% imSurf.x_s (x-axis vector in pixels), imSurf.z_s (surface elevation vector in pixels
% (upside down, upright surface is 2048 - imSurf.z_s + 1))
%% Author:
% Marc Buckley
%% Last update:
% 06/24/2013
%% Example:
% load('E:\data\Exp4\RawImages\Lfv\Exp4_Lfv_1004.mat')
% imSurf = findSurface(img);
% figure, imagesc(img), colormap(bone), hold on, plot(imSurf.x_s, imSurf.z_s, 'r')
%%
function imSurf = findSurface_lc(img)

surface = nan(1,size(img,2));
% 
for i=1:size(img,2)
    imgi = img(:,i);
    %locate outliers and nan them
    imgistd = std(imgi);
    imgimean = mean(imgi);
    imgi(abs(imgi-imgimean)>3*imgistd)=nan;
    %smooth each column
    smth_imgi = smoothn(imgi,10000);
    %compute gradient on each column
    gv =  gradient(smth_imgi);
    %find first large peak in gradient
    [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/3, 'npeaks',1);
    surface(i) = locs;
    clear imgi gv locs pks
end
%
surface=despike_fab(surface);
%locate outliers in surface and nan them
surfstd = std(surface);
surfmean = mean(surface);
surface(abs(surface-surfmean)>3*surfstd)=nan;
smth_surface= smoothn(surface,1000);

imSurf.x_s = 1:size(img,2);
imSurf.z_s = smth_surface;
imSurf.z_s_raw = surface;
imSurf.img=img;



%%
% imgi=img;
% imgistd=std2(img);
% imgimean=mean2(img);
% imgi(abs(imgi-imgimean)>3*imgistd)=nan;
% smth_imgi = smoothn(imgi,10000);
% [gvx,gvy]=gradient(smth_imgi);
% for i=1:size(img,2)
% [~,locs] = findpeaks(gvy(:,i),'minpeakheight', max(gvy(:,i))/2, 'npeaks',1);
% surface(i) = locs;
% end

