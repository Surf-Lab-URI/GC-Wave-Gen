clear
clc
close all

% Define Path
DataPath = '/media/surflab/New Volume/ExpPilot/ExpPilot5/ExpPilot5_Scene2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_Andy/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water
%%
frames = 0:1:7049; %Range of frame numbers (the numbers in the file names, which start at zero).
nF = length(frames);
fXs = zeros(nF, 4106);
fYs = zeros(nF, 4106);

pps = 14.5; %pairs per second

spp = 1/pps; % seconds per pair

yLimits = [1950,2150];
dy = yLimits(2)-yLimits(1);

dt_pair = 22.222e-3; % time between pictures in a given pair

dydt = dy/dt_pair;

DeltaT = nF/2*spp;


surfSigmas = [50,40,30,20,10];
surfSteps = [50,40,30,20];
surfMask = 1;
%% Assemble array of times since 
t = zeros(1,nF);
t(1) = floor(frames(1)/2)*spp + mod(frames(1),2)*dt_pair;
for i = 2:nF
    if mod(frames(i),2)==1
        t(i) = t(i-1)+dt_pair;
    else
        t(i) = t(i-1)+spp-dt_pair;
    end
end
%%
tic
parfor i = 1:nF
    idx = frames(i);

    PIV1Dir_temp = PIVWaterDir;

    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    % imagename
    
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';

    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

    % imagesc(PIVSurfW1_CamAngle, [0,70])
    % hold on
    % set(gca,'DataAspectRatio',[1 1 1])
    % ylim([1900,2200])

    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface] = FindWaterSurface(PIVSurfW1_CamAngle, surfSigmas, surfSteps, surfMask);

    
    % plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
    % pause
    
    fXs(i,:) = XPIVSurfW1_Surface;
    fYs(i,:) = PIVSurfW1_Surface;
end
toc
%%
figure(2)
for i = 1:nF
    plot(fXs(i,:), -fYs(i,:))
    set(gca,'DataAspectRatio',[1 1 1])
    ylim([-2200,-1900])
    pause(0.5)
end

%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle,surfSigmas, surfSteps, surfMask)

    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = CrapperOptimized_FindSurface(PIVSurfW_CamAngle(1:2800,X), surfSigmas, surfSteps, surfMask);
    PIVSurf_Surface_Raw = imSurf.surface;

    % f = gcf;
    % figure(f)
    % hold on
    % plot(X,PIVSurf_Surface_Raw,'-g')
    
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);

    % plot(X,PIVSurf_Surface_Int,'-k')


    PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
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

