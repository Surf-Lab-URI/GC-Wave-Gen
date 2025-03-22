clear
close all
clc

ROOTPath = '/Volumes/New Volume/ExpPilot/'; % FabioASI

ExpPilots = [1];%[ 1 1 2 2 3 3 4 4 4 5 5 ];
Scenes = [2];%[ 1 2 1 2 1 2 1 2 3 1 2 ] ;

frames = 460;

nF = length(frames);
fXs = zeros(nF, 8050);
fYs = zeros(nF, 8050);

mpp = 6.2615e-05; %  (meters per pixel for the LFV air camera)

for i = 1:length(ExpPilots) % Main Loop
    
    expName = num2str(ExpPilots(i));
    
    sceneName = num2str(Scenes(i));
    
    % Define Path
    DataPath = [ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
    LoadPath = [DataPath 'RAW/'];
    RawDataPath = [DataPath 'RAW/'];
    ResultsPath = [DataPath 'RESULTS_Andy/'];
    
    % Air
    AirPath = [ResultsPath 'Air/'];
    SavePIVAirPath = [AirPath 'PIV_Velocities_raw/'];
    SaveSurfAirPath = [AirPath 'Surfaces/'];
    FieldsAirPath  = [AirPath 'CALCULATED_FIELDS/'];
    
    % Water
    WaterPath = [ResultsPath 'Water/'];
    SavePIVWaterPath = [WaterPath 'PIV_Velocities_raw/'];
    SaveSurfWaterPath = [WaterPath 'Surfaces/'];
    FieldsWaterPath  = [WaterPath 'CALCULATED_FIELDS/'];
    
    if ~exist(SavePIVAirPath, 'dir')
        mkdir(SavePIVAirPath);
    end
    if ~exist(SavePIVWaterPath, 'dir')
        mkdir(SavePIVWaterPath);
    end
    
    if ~exist(SaveSurfAirPath, 'dir')
        mkdir(SaveSurfAirPath);
    end
    if ~exist(SaveSurfWaterPath, 'dir')
        mkdir(SaveSurfWaterPath);
    end
    
    %% Bad Frame, Normalization and lens distortion parameters 
    load cameraParams.mat
    
    %% Parameters    
    IntrWndw_A = [128 64 32 16 8]; %[64 32 16 8]; %Interogation Window (size of box of pixels used for cross-correlation)
    GrdSpc_A = [64 32 16 8 4]; %[32 16 8 4];
    IntrWndw_W = [[256 64 32 16 8]*8 16 8];
    GrdSpc_W = [[128 32 16 8 4]*8 8 4];
    
    CST.DX = 38.106d-06; % meters per pixel (~38 micron/pix) in PIV Air pixel resolution
    CST.DX_W = 41.782d-6; % meters per pixel (~38 micron/pix) in PIV Water pixel resolution
%     CST.DT_A = DeltaT_A(i); %Delta t for air (s)
    CST.DT_W = 22.22222d-3; %Delta t for water (s)
    CST.GS = GrdSpc_A(end);
    CST.IW = IntrWndw_A(end);
    
    
    %%% Save Parameters
    save([ResultsPath 'ExpPilot' expName '_Scene' sceneName '_Parameters.mat'],'CST')
    
    %% Frame to process
    PIVAirDir = dir([LoadPath 'PIV Air/' '*.raw']); %Find all the raw air files and store their names and directory in a struct
    PIVWaterDir = dir([LoadPath 'PIV Water/' '*.raw']); %Same for water
    FI = 0;%First index
    LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1); %Last index. Stop one short of last frame bc they are processed in pairs.
    
    %% Processing frames
    image_index = FI+1:2:LI; %1, 3, 5,... Set of indices to loop through. Images are processed in pairs, hence the increment of 2
    framesCtr = 1;
    for idx = frames%:numel(image_index) % Main Loop through the pairs of images. Starting on the 10th pair of images for some reason.
        
        % Indexes for images
        pair_index = (image_index(idx)+1)/2; %number this pair of images
%         PIV1Dir_temp = PIVAirDir;
        SurfDir_temp = PIVAirDir;
%         ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4); %number of the first image in the pair as it appears in the raw file name (i.e. starting from 0000)
%         ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);%number of the second image in the pair as it appears in the raw file name
%         
%         ImageNum_Water1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4); %number of the first image in the pair as it appears in the raw file name (i.e. starting from 0000)
%         ImageNum_Water2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);%number of the second image in the pair as it appears in the raw file name
        PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4); %number of the pair as it appears in the file name
        
%         if i == 1 % I don't understand this
%             ImageNum_Air1 = PIV1Dir_temp(image_index(idx)-1).name(max(strfind(PIV1Dir_temp(image_index(idx)-1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
%             ImageNum_Air2 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
%         end
        
        %% Air
        
              
        %%% Surface detection
        % Load PIVSurf Air - LFV
        imagename = [LoadPath 'PIVSurf Air - LFV/ExpPilot' expName '_Scene' sceneName' '_PIVSurf Air - LFV_' PairNum '.raw'];
        [IM1] = (load_Image_IOCoreView_48MP(imagename));
        PIVSurf_A_Raw = IM1;
        PIVSURF_A = PIVSurf_A_Raw./(smooth(mean(PIVSurf_A_Raw(2400:end,:)),1000)/max(smooth(mean(PIVSurf_A_Raw(2400:end,:)),1000)))'; %smoothing and scaling
        
        % Lens distortion correction
        [PIVSurfA_Undistorted] = PIVSurfA_LFV_LensDistCorr(IM1,cameraParams);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE DONE!!!!!
        
        % Camera Angle Correction
        [PIVSurfA_CamAngle] = PIVSurfAir_LFV_CamAngle_Correction(PIVSurfA_Undistorted);
        
        % Extract surface
        [BadFrame, XLFV_Surface, LFV_Surface] = FindSurface(PIVSurfA_CamAngle);
        
        figure(1)
        plot(mpp*XLFV_Surface, -mpp*LFV_Surface)
        xlabel('x (m)')
        ylabel('y (m)')
        axis equal
        pause(0.01)
        
        fXs(framesCtr,:) = XLFV_Surface;
        fYs(framesCtr,:) = LFV_Surface;
        
        framesCtr = framesCtr + 1;
    end
end
save('surfaceTestAndy.mat','fXs','fYs');

%%
is = length(fXs(:,1));
hold off
for i = 100:is
    figure(1)
    plot(mpp*fXs(i,:), -mpp*fYs(i,:))
    hold on
    xlabel('x (m)')
    ylabel('y (m)')
    set(gca,'DataAspectRatio',[1 1 1])
    ylim([-0.2,-0.1])
    xlim([0,0.6])
    pause
    hold off
end

%%
%% Step through raw frames
clear
clc


frames = 1090:1:1120;

nF = length(frames);
fXs = zeros(nF, 4106);
fYs = zeros(nF, 4106);

i = 1;


% Define Path
DataPath = '/Volumes/New Volume/ExpPilot/ExpPilot5/ExpPilot5_Scene2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_Andy/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water


FI = 0;%First index
LI = length(PIVWaterDir)-1;
image_index = FI+1:LI; %1, 3, 5,... Set of indices to loop through. Images are processed in pairs, hence the increment of 2
framesCtr = 1;

for idx = frames%:numel(image_index) % Main Loop through the pairs of images. Starting on the 10th pair of images for some reason.
    idx
    PIV1Dir_temp = PIVWaterDir;
%     ImageNum_Water1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4); %number of the first image in the pair as it appears in the raw file name (i.e. starting from 0000)

%     imagename = [LoadPath 'PIVSurf Water/ExpPilot' expName '_Scene' sceneName' '_PIVSurf Water_' ImageNum_Water1 '.raw'];
    
    imagename = [PIV1Dir_temp(idx).folder '/' PIV1Dir_temp(idx).name];
    
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';

%     [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
%     [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

%     [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface] = FindWaterSurface(PIVSurfW1_CamAngle);
    
    figure(20)
    imagesc(PIVSURF_W1/mean(PIVSURF_W1,'all'), [0,3])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    ylim([1900,2200])
%     plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
    pause
    
%     fXs(framesCtr,:) = XPIVSurfW1_Surface;
%     fYs(framesCtr,:) = PIVSurfW1_Surface;

    framesCtr = framesCtr + 1;
end

%%
is = length(fXs(:,1));
hold off
writerObj = VideoWriter('ExpPilot5Scene2Clip');
writerObj.FrameRate = 5;
open(writerObj);
for i = 1:is
    figure(21)
    plot(fXs(i,:), -fYs(i,:)) 
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    ylim([-2200,-1900])
    title(LoadPath)
    hold off
    
    frame = getframe(gcf)   ; 
    writeVideo(writerObj, frame);    
    pause(0.2) 
end

close(writerObj)

%%
function [BadFramePIVSurfLFV,XLFV_Surface,LFV_Surface] = FindSurface(PIVSurfA_CamAngle)
%% Step 1: Find surface
% Extract surface
XX = 1001:3000;
YY = 251:8300;
% [imSurf] = Copy_of_FindSurface(PIVSurfA_CamAngle(1001:3000,:), 5, 5);
[imSurf] = Copy_of_FindSurface(PIVSurfA_CamAngle(XX,YY), 5, 5);
PIVSurf_Surface_Raw = imSurf.surface;
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1d4); %10000); %Generates a spline
if length(SP.coefs)>2
    PIVSurf_Surface_A = SP.coefs(2:end-1);
else
    A = polyfit([YY(1) YY(end)], [SP.coefs(1) SP.coefs(end)],1);
    PIVSurf_Surface_A = polyval(A,YY);
end
% Usurf = 501:length(PIVSurf_Surface_A)-500;
% Vsurf = PIVSurf_Surface_A(501:end-500)+1000;
Usurf = YY;
Vsurf = PIVSurf_Surface_A+XX(1)-1;

%% Check if bad frame
BadFramePIVSurfLFV = 0;
if imSurf.badFrameBool == 1
    BadFramePIVSurfLFV = 1;
end
%%
XLFV_Surface = Usurf;
LFV_Surface = Vsurf;

end
%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle)

%% PIVSurf Water surface detection

%% Step 1: Find surface
% Extract surface
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