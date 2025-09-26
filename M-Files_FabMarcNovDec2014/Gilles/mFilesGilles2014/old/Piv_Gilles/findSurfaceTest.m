% tic
% clear all
% close all
% load('D:\Gilles_LC_PIV\surfaceDetection\Exp4_Lfv_0050.mat')
% img=imgLfv;
% mask = nan(size(img));
% surface = nan(1,size(img,2));
% 
% for i=1:size(img,2)
%     j=1;
%     while img(j,i)<1000
%         j=j+1;
%     end
%     surface(i)=j;
% end
% surface=smoothn(surface,10);
% figure,imagesc(img), colormap(bone)
% hold on, plot(surface,'r')
% toc
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imSurf = findSurfaceLfv(imgLfv,3);   %!!! Marche mieux avec 3 9 (init 2)
% hold on, plot(imSurf.x_s, imSurf.z_s, 'r')
% toc
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear all
close all
load('D:\Gilles_LC_PIV\surfaceDetection\Exp4_Lfv_0000.mat')
img=imgLfv;
surfEdge=edge(imgLfv,'sobel');
surfEdgeSmt=5*smoothn(surfEdge);
for i=1:size(surfEdgeSmt,2)
    j=1;
    while surfEdgeSmt(j,i)<0.9
        j=j+1;
    end
    surface(i)=j;
end
surface=smoothn(surface,100);
figure,imagesc(img), colormap(bone)
hold on, plot(surface,'g')
toc






tic
clear all
close all
load('D:\Gilles_LC_PIV\surfaceDetection\Exp4_Lfv_0000.mat')
img=imgLfv;
imSurf = findSurfaceLfv(imgLfv,3);   %!!! Marche mieux avec 3 (init 2)
figure,imagesc(img), colormap(bone)
hold on, plot(imSurf.x_s, imSurf.z_s, 'r')
toc


% surfEdge=edge(imgLfv,'prewitt');
% 'prewitt' 'sobel'  



% 
% 
% imSurf = findSurfaceLfv(surfEdgeSmt,3);
% figure,imagesc(img), colormap(bone)
% hold on , imshow(surfEdgeSmt),colorbar
% hold on, plot(imSurf.x_s, imSurf.z_s, 'r')
% toc
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 

















