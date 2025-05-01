clear
clc

%%
DataPath = '/media/surflab/New Volume/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water


PIV1Dir_temp = PIVWaterDir;

mpp = 6.236119402985075e-5; %Accurate only for the main dataset, not the pilot

% -ExpAW1Run2 around 1340 but not as prominant
% -ExpAW5 Run2 super steep parasitic? capillaries. Perhaps periodic crapper
% capillaries. solitons? starting around 824. Steep soliton around 861.
% Some good asymmetry

%% Flip through frames with minimal processing and scale bar
viewing = true;

load Norm_PIV.mat

idx = 0;
while viewing

    idx
    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    [IM1] = load_Image_IOCoreView_12MP(imagename);

    % if mod(idx, 2) == 0
    %     IM1 = IM1./Norm_PIVSurfW1;
    % else
    %     IM1 = IM1./Norm_PIVSurfW2;
    % end

    IM1 = IM1/mean(IM1(:),'omitnan')*20;
    % IM1(IM1>1023) = 1023;
    
    figure(50)
    hold off
    imagesc(IM1,[0,70])
    hold on
    colormap gray
    set(gca,'DataAspectRatio',[1 1 1])
    axis off
    set(gca,'FontSize',24)
    % title(s,'Interpreter','none')
    
    % ylim([1750,2150])
    
    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/mpp;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.05)*[1 1];
    try
        delete(sb)
        delete(sbt)
    end
    sb = plot(xsb,ysb,'-k', 'LineWidth',10);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex');    
    
    % delete(p2)

    ip = input("Next Frame?","s")
    ip

    if ip == 'a'
        idx = max(idx-1,0);
    elseif ip == 'd'
        idx = min([size(PIVWaterDir,1)-1,idx+1]);
    else
        try
            ip = uint16(str2double(ip));
            idx = ip;
        end
    end
end
%% Flip through frames with edge detection
load cameraParams.mat
load Norm_PIV.mat
load BadPix.mat

viewing = true;

load Norm_PIV.mat

idx = 0;
while viewing

    idx
    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    [IM1] = load_Image_IOCoreView_12MP(imagename);

    % if mod(idx, 2) == 0
    %     IM1 = IM1./Norm_PIVSurfW1;
    % else
    %     IM1 = IM1./Norm_PIVSurfW2;
    % end

    IM1 = IM1/mean(IM1(:),'omitnan')*15;

    for iiii = 1:length(BadPix_SurfW)
        IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
    end
    
    PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
    % PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;
    
    
    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSurf_W1_Raw);
    
    % Camera Angle Correction
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);
    
    PIV1_W = [];
    % Extract surface
    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');

    
    figure(50)
    hold off
    imagesc(PIVSurfW1_CamAngle,[0,70])
    hold on
    colormap gray
    set(gca,'DataAspectRatio',[1 1 1])
    axis off
    set(gca,'FontSize',24)

    plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'-r','LineWidth',2)
    xlim([XPIVSurfW1_Surface(1),XPIVSurfW1_Surface(end)])
    ylim([1750,2250])
    
    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/mpp;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb)
        delete(sbt)
    end
    sb = plot(xsb,ysb,'-k', 'LineWidth',10);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex');    
    
    % delete(p2)

    ip = input("Next Frame?","s")
    ip

    if ip == 'a'
        idx = max(idx-1,0);
    elseif ip == 'd'
        idx = min([size(PIVWaterDir,1)-1,idx+1]);
    else
        try
            ip = uint16(str2double(ip));
            idx = ip;
        end
    end
end

%% Save frames for a video
clear f F

load cameraParams.mat
load Norm_PIV.mat
load BadPix.mat

viewing = true;

load Norm_PIV.mat

idx = 0;

figure('units','pixels','Position',[0,0,1500,200])

f = 1;
tic
idxs = 800:2:900
F = struct('cdata',cell(length(idxs),1),'colormap',cell(length(idxs),1));
parfor i = 1:length(idxs)
    idx = idxs(i)
    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    [IM1] = load_Image_IOCoreView_12MP(imagename);

    % if mod(idx, 2) == 0
    %     IM1 = IM1./Norm_PIVSurfW1;
    % else
    %     IM1 = IM1./Norm_PIVSurfW2;
    % end

    IM1 = IM1/mean(IM1(:),'omitnan')*15;

    for iiii = 1:length(BadPix_SurfW)
        IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
    end
    
    PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
    % PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;
    
    
    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSurf_W1_Raw);
    
    % Camera Angle Correction
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);
    
    PIV1_W = [];
    % Extract surface
    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');

    hold off
    imagesc(PIVSurfW1_CamAngle,[0,70])
    hold on
    colormap gray
    set(gca,'DataAspectRatio',[1 1 1])
    axis off
    set(gca,'FontSize',24)

    plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'-r','LineWidth',2)
    xlim([XPIVSurfW1_Surface(1),XPIVSurfW1_Surface(end)])
    ylim([1750,2250])
    
    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/mpp;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb)
        delete(sbt)
    end
    sb = plot(xsb,ysb,'-k', 'LineWidth',10);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex');
    F(i) = getframe(gcf);

    fname = ['videoframes/test' num2str(i)]; % full name of image
    print('-djpeg','-r600',fname)     % save image with '-r200' resolution
end
toc
%% Generate Video File from frames save in previous section
vw = VideoWriter('test.avi', 'Uncompressed AVI');
vw.FrameRate = 5;
open(vw);
vw

for i = 1:length(idxs)
    fname = ['videoframes/test' num2str(i)]; % full name of image
    I = imread([fname '.jpg']);       % read saved image
    frame = im2frame(I);              % convert image to frame
    writeVideo(vw,frame)
end

close(vw);

%% Test Surface Detection and transformation
load cameraParams.mat
load Norm_PIV.mat
load BadPix.mat

idx = 862;
imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];

[IM1] = load_Image_IOCoreView_12MP(imagename);

% figure(50)
% imagesc(IM1,[0,70])
% colormap gray
% set(gca,'DataAspectRatio',[1 1 1])

% Remove Bad Pixels and interpolate
for iiii = 1:length(BadPix_SurfW)
    IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
end
PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;


[PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % TO BE DONE!!!!! BUT PROBABLY NOT NEEDED

                % Camera Angle Correction
[PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);


PIV1_W = [];
                % Extract surface
[BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');
