%% Name of the experiments

ROOTPath = '\\spray3\d\data\EXPERIMENTS\'; % FabioASI

ExpAW = 8;
Run = 1;
Acc = 0.13;
Wind = 7;
DeltaT_A = 90d-6;
i = 1;

%% Folders and parameters
expName = ['ExpAW' num2str(ExpAW(i)) '_acc' num2str(Acc(i)) '_W' num2str(Wind(i)) 'V' ];

if ExpAW(i) == 6
    expName = ['ExpAW' num2str(ExpAW(i)) '_acc' num2str(Acc(i)) '_W' num2str(Wind(i)) 'V_LidOpen' ];
end

runName = [ 'Run' num2str(Run(i))];

% Define Path
DataPath = [ROOTPath expName '\' expName '_' runName '\' ];
LoadPath = [DataPath 'RAW\'];
RawDataPath = [DataPath 'RAW\'];
ResultsPath = [DataPath 'RESULTS_fabio\'];

% Air
AirPath = [ResultsPath 'Air\'];
SavePIVAirPath = [AirPath 'PIV_Velocities_raw\'];
SaveSurfAirPath = [AirPath 'Surfaces\'];
FieldsAirPath  = [AirPath 'CALCULATED_FIELDS\'];

% Water
WaterPath = [ResultsPath 'Water\'];
SavePIVWaterPath = [WaterPath 'PIV_Velocities_raw\'];
SaveSurfWaterPath = [WaterPath 'Surfaces\'];
FieldsWaterPath  = [WaterPath 'CALCULATED_FIELDS\'];

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

%%% Bad Frame, Normalization and lens distortion parameters
load cameraParams.mat
load Norm_PIV.mat

%%% Parameters
IntrWndw_A = [128 64 32 16 8]; %[64 32 16 8];
GrdSpc_A = [64 32 16 8 4]; %[32 16 8 4];
IntrWndw_W = [[256 192 128 96 64 48 32 16 8]*8 32 16 8]; %[[128 64 32 16 8]*8 16 8];
GrdSpc_W = IntrWndw_W/2; %[[64 32 16 8 4]*8 8 4];

CST.DX = 40.126d-06; % meters per pixel (~40 micron/pix) in PIV Air pixel resolution
CST.DX_W = 41.242d-6; % meters per pixel (~41 micron/pix) in PIV Water pixel resolution
CST.DT = DeltaT_A(i);
CST.DT_W = 22.22222d-3;
CST.GS = GrdSpc_A(end);
CST.IW = IntrWndw_A(end);

%%% Frame to process
PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
FI = 0;
LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
image_index = FI+1:2:LI;

%% Water analysis
count = 40;
for idx = 10:101:1020
    count = count+1;
    pair_index = (image_index(idx)+1)/2;
    PIV1Dir_temp = PIVAirDir;
    SurfDir_temp = PIVAirDir;
    ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
    ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);

    ImageNum_Water1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
    ImageNum_Water2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);

    PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4);
    % First PIV image Water
    filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water1 '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM1_W = IM1;

    % Second PIV image Water
    filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water2 '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM2_W = IM1;

    % Camera Angle Correction
    [PIVWater1_CamAngle] = PIVWater_CamAngle_Correction(IM1_W); %Image 1
    [PIVWater2_CamAngle] = PIVWater_CamAngle_Correction(IM2_W); %Image 2
    PIVWater1_CamAngle(PIVWater1_CamAngle<0) = 0;
    PIVWater2_CamAngle(PIVWater2_CamAngle<0) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-process
    %%% PIVWater 1
    PIV = PIVWater1_CamAngle;
    [PIV] = Pre_process_PIV_Image_Water_IM1(PIV); %pre_proc 1st image
    PIV1_W = PIV;
    %%% PIVWater 2
    PIV = PIVWater2_CamAngle;
    [PIV] = Pre_process_PIV_Image_Water_IM2(PIV); %pre_proc 2nd image
    PIV2_W = PIV;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface detection
    %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W1
    imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water1 '.raw'];
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW;

    % Lens distortion correction
    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE
    % DONE!!!!! BUT PROBABLY NOT NEEDED

    % Camera Angle Correction
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

    % Extract surface
    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');

    % Mask PIVWater1
    [Mask1_W] = PIVWater_Mask(PIV1_W, PIVW1_Surface);

    %%%%%% COORDINATE TRANSFORMATION DATA
    [h, w] = size(PIV1_W);
    [Uorb_W1,Vorb_W1] = GenerateTransfo_Fabio_airwater(XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, CST, h, w, PIV1_W);

    %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W2
    imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water2 '.raw'];
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W2_Raw = IM1;
    PIVSURF_W2 = PIVSurf_W2_Raw./Norm_PIVSurfW;

    % Lens distortion correction
    [PIVSurfW2_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W2);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE
    % DONE!!!!! PROBABLY NOT NEEDED

    % Camera Angle Correction
    [PIVSurfW2_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW2_Undistorted);

    % Extract surface
    [BadFramePIVSurfW2, XPIVSurfW2_Surface, PIVSurfW2_Surface, XPIVW_PIVSurfW2_Surface, PIVW_PIVSurfW2_Surface, PIVW2_Surface,T2_W,RotAngle_W,DY_W] = ExtractSurface_PIVSurfWater(PIVSurfW2_CamAngle,PIV2_W,'2');

    %         % Transform surface to mask PIVWater Surface
    %         [XPIVWater2_Surface,PIVWater2_Surface,XPIV_PIVWater2_Surface_CRR,PIV_PIVWater2_Surface_CRR] = Correct_Surface_Water(PIVWater2_CRR,XPIVW_PIVSurfW2_Surface,PIVW_PIVSurfW2_Surface,PIV2_W,idx);

    % Mask PIVWater2
    [Mask2_W] = PIVWater_Mask(PIV2_W, PIVW2_Surface);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIV CALCULATIONS
    [CompVelWater] =  ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIV1_W, PIV2_W, Mask1_W, Mask2_W, Mask1_W, IntrWndw_W, GrdSpc_W, Uorb_W1,Vorb_W1);
    figure(count);imagesc(CompVelWater.INTdelx.*CompVelWater.Mask);title('INTdelx')
    figure(count+100);imagesc(CompVelWater.dcor);clim([0 1]);title('dcor')
end