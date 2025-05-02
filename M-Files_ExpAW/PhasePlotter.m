%% Step through raw frames
clear
clc
close all

% Define Path
DataPath = '/media/surflab/New Volume/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_Andy/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water
%%
frames = 850:1:870;
nF = length(frames);
fXs = zeros(nF, 3639);
fYs = zeros(nF, 3639);

pps = 14.5; %pairs per second

spp = 1/pps; % seconds per pair

yLimits = [1950,2150];
dy = yLimits(2)-yLimits(1);

dt_pair = 22.222e-3; % time between pictures in a given pair

dydt = dy/dt_pair;

DeltaT = nF/2*spp;

% FI = 0;%First index
% LI = length(PIVWaterDir)-1;
% image_index = FI+1:LI; %1, 3, 5,... Set of indices to loop through. Images are processed in pairs, hence the increment of 2

%%
CompImg = NaN(ceil(DeltaT*dydt),4176);
framesCtr = 1;
for idx = frames
    idx
    PIV1Dir_temp = PIVWaterDir;

    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    imagename
    
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';

    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

    PIV1_W = [];
    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');
    
    if mod(idx,2) == 0
        CompImg(round((framesCtr-1)/2*spp*dydt)+1:round((framesCtr-1)/2*spp*dydt)+dy+1,1:size(PIVSurfW1_CamAngle,2)) = PIVSurfW1_CamAngle(yLimits(1):yLimits(2),:)/mean(PIVSurfW1_CamAngle,'all');
    else
        CompImg(round(((framesCtr-2)/2*spp+dt_pair)*dydt)+1:round(((framesCtr-2)/2*spp+dt_pair)*dydt)+dy+1,1:size(PIVSurfW1_CamAngle,2)) = PIVSurfW1_CamAngle(yLimits(1):yLimits(2),:)/mean(PIVSurfW1_CamAngle,'all');
    end

    imagesc(PIVSurfW1_CamAngle, [0,70])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    ylim([1900,2200])
    plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
    pause(0.1)
    
    fXs(framesCtr,:) = XPIVSurfW1_Surface;
    fYs(framesCtr,:) = PIVSurfW1_Surface;

    framesCtr = framesCtr + 1;
end
%%
P_threshold = 10;

figure(2)

mpp = 6.493178e-5; %for main dataset only
x = (1:size(CompImg,2))*mpp;
y = 1:size(CompImg,1);
t = y/dydt;

imagesc(CompImg,'XData',x,'YData',t,[0,3])
hold on
xlabel('x (m)')
ylabel('t (s)')
s = DataPath(end-23:end-1) + " " + frames(1) + " to " + frames(end);
title(s,'Interpreter','none')
set(gca,'DataAspectRatio',[1*mpp 1/dydt 1])
set(gca, 'ytick', 0:spp:t(end));
set(gca, 'xtick', 0:0.02:x(end));
set(gca,'FontSize',24)
colormap gray
% axis off

minPts = NaN(300,2);
cmp = 1;
for i = 1:nF
    [TF, P] = islocalmin(-fYs(i,:));
    minsMask = P > P_threshold;

    if mod(i,2) == 1
        offset = round((i-1)/2*spp*dydt)+1-yLimits(1);
        offset2 = round((i-1)/2*spp*dydt)+1;
        plot(fXs(i,:)*mpp,(fYs(i,:)+offset)/dydt,'r')
        plot(fXs(i,minsMask)*mpp,(fYs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
    else
        offset = round(((i-2)/2*spp+dt_pair)*dydt)+1-yLimits(1);
        offset2 = round(((i-2)/2*spp+dt_pair)*dydt)+1;
        plot(fXs(i,:)*mpp,(fYs(i,:)+offset)/dydt,'r')
        plot(fXs(i,minsMask)*mpp,(fYs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
    end
    
    minPts(cmp:cmp+length(fXs(i,minsMask))-1,:) = [fXs(i,minsMask)',offset2*ones(1,length(fXs(i,minsMask)))'];
    cmp = cmp + length(fXs(i,minsMask));
    hold on
    % pause
end
minPts(any(isnan(minPts), 2), :) = [];


% plot(minPts(:,1),minPts(:,2),'r*') 

%%
troughImg = zeros(size(CompImg));
idxs = sub2ind(size(troughImg), minPts(:,2),minPts(:,1));
troughImg(idxs) = 1;
troughImg = imgaussfilt(troughImg,10);
figure(5)

troughImg = (troughImg > 0.00001);
imagesc(troughImg)

figure(3)
[H,theta,rho]= hough(troughImg,'Theta',-45:0.01:-10);
% Hblur = imgaussfilt(H,5);
imshow(imadjust(rescale(H)),'XData',theta,'YData',rho,...
      'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);

P_threshold = 0.15;
[TF1, P1] = islocalmax(Hblur,1);
[TF2, P2] = islocalmax(Hblur,2);
P = P1.*P2;
figure(4)
imagesc(P,'XData',theta,'YData',rho)
[iPr,iPt] = find(P > P_threshold);
thetas = theta(iPt);
rhos = rho(iPr);

figure(1)
hold on
ms = NaN(1,length(thetas));
for i = 1:length(thetas)
    x = get(gca,'XLim');
    ms(i) = -cos(theta(i))/sin(theta(i));
    y = (rhos(i) - x* cos(theta(i)))/ sin(theta(i));
    plot(x,y,'r');
end
% plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
%%
% plotLine([1734,1366], [1978,1571],1)
plotLine([2091,1972],[2335,2169],1)
plotLine([242,1963],[477,2163],1)
plotLine([1621,923],[1388,725],1)
plotLine([327,1341],[574,1544],1)

%%
function plotLine(A,B,figNum)
    figure(figNum)
    xlim = get(gca,'XLim');
    m = (B(2)-A(2))/(B(1)-A(1));
    b = B(2) - m*B(1);
    y1 = m*xlim(1) + b;
    y2 = m*xlim(2) + b;
    hold on
    plot([xlim(1) xlim(2)],[y1 y2],'-r','LineWidth',2)
    plot([A(1) B(1)],[A(2) B(2)],'*r','MarkerSize',20)
    hold off
end
%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle)
    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = Copy_of_FindSurface(PIVSurfW_CamAngle(1:2800,X), 5, 5);
    PIVSurf_Surface_Raw = imSurf.surface;
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
    [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1d3);
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