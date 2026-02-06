clear
clc
close all

% This script is for finding the upwind and downwind ends of the PIV
% image in the IR image. The results I got using this script are given
% below. These points are coordinates in the transformed IR images:
%           x           y
% Upwind    291.58      202
% Downwind  289.13      511

%% Load image showing the location of the PIV laser sheet in the IR frame
calimgname = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/IRcalib/calib/PIV laser LON.png';

a=imread(calimgname);
size(a)
figure
imagesc(a)
a = a(1:512,1:640,:);
a = rgb2gray(a);
size(a)
figure(1)
imagesc(a)

%%
a=medfilt2(double(a));
U = [94 25; 480 22; 51 500; 495 501];
X = [51 22; 495 22; 51 501; 495 501];
T = maketform('projective', U, X); 
[b1, ~, ~] = imtransform(a,T,'XYScale',1); % Rectification for IR
b1=b1(7:518,64:698);  %croping
imgIR=imresize(b1,[635 635]); % making calibration same in x and z

PIV_laser_LON = imgIR;
IR.img=imgIR;
IR.DX=3.87e-004;

figure(2)
imagesc(IR.img)

%% Select points at the corners of the long skinny rectangle that indicates the position of the laser sheet
figure(2)
imagesc(IR.img)
disp('click on top corners')
[xtop,ytop] = ginput(2);
xtop
ytop
disp('click on bottom corners')
[xbot,ybot] = ginput(2);
xbot
ybot
%%
mean(xtop)
mean(xbot)
diff(xtop)*IR.DX

%% Find downwind edge of PIV image
dwcalimgname = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/IRcalib/calib/downwind edge of PIV (in laser sheet).png';
a=imread(dwcalimgname);
size(a)
figure(3)
hold off
% imagesc(a)
a = a(1:512,1:640,:);
a = rgb2gray(a);
size(a)
a=medfilt2(double(a));
T = maketform('projective', U, X); 
[b1, ~, ~] = imtransform(a,T,'XYScale',1); % Rectification for IR
b1=b1(7:518,64:698);  %croping
imgIR=imresize(b1,[635 635]); % making calibration same in x and z

dwImg = imgIR;
IR.img=imgIR;
IR.DX=3.87e-004;

figure(3)
imagesc(IR.img)
%choose circle coordinates and radius that match the intersection of the
%bob line with the water
hold on
cx = 243; cy = 520; r = 9;
circle(cx,cy,r)
dwPIVbob = [cx, cy - r];
plot(dwPIVbob(1),dwPIVbob(2),'*','MarkerSize',10)

%% Find upwind edge of PIV image
uwcalimgname = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/IRcalib/calib/upwind edge of PIV (in laser sheet).png';
a=imread(uwcalimgname);
size(a)
figure(4)
hold off
% imagesc(a)
a = a(1:512,1:640,:);
a = rgb2gray(a);
size(a)
a=medfilt2(double(a));
T = maketform('projective', U, X); 
[b1, ~, ~] = imtransform(a,T,'XYScale',1); % Rectification for IR
b1=b1(7:518,64:698);  %croping
imgIR=imresize(b1,[635 635]); % making calibration same in x and z

uwImg = imgIR;
IR.img=imgIR;
IR.DX=3.87e-004;


imagesc(IR.img)
hold on
uwPIVbob = [278,202];
plot(uwPIVbob(1),uwPIVbob(2),'*','MarkerSize',10)
%% Determine upwind and downwind edges of the PIV image
ybm = mean(ybot);
ytm = mean(ytop);
xbm = mean(xbot);
xtm = mean(xtop);
m = (ytm-ybm)/(xtm-xbm);
uwPIV = [1/m*(uwPIVbob(2)-ybm) + xbm, uwPIVbob(2)];
dwPIV = [1/m*(dwPIVbob(2)-ybm) + xbm, dwPIVbob(2)];
%% plot upwind and downwind edges and line all together
figure(5)
imagesc(PIV_laser_LON)
hold on
plot([mean(xtop),mean(xbot)],[mean(ytop),mean(ybot)],'-k')
plot(uwPIVbob(1),uwPIVbob(2),'*r','MarkerSize',10)
plot(dwPIVbob(1),dwPIVbob(2),'*k','MarkerSize',10)
plot(uwPIV(1),uwPIV(2),'.r','MarkerSize',10)
plot(dwPIV(1),dwPIV(2),'.k','MarkerSize',10)
%% Function for drawing circles
function h = circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
end