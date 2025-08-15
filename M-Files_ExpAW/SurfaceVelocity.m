clear
clc

ROOTPath = '/media/surflab/Working24/ExpAW/';

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

i = 5; % ExpNumber
runNum = 2;

ExpAW = ExpDir(i).name(6);
Acc = ExpDir(i).name(11:14);
Wind = ExpDir(i).name(17);
switch Wind
    case '4'
        DeltaT_A = 200d-6;
    case '5'
        DeltaT_A = 200d-6;
    case '6'
        DeltaT_A = 120d-6;
    case '7'
        DeltaT_A = 90d-6;
end


expName = ['ExpAW' ExpAW '_acc' Acc '_W' Wind 'V' ];

if strcmp(ExpAW,'6')
    expName = ['ExpAW' ExpAW '_acc' Acc '_W' Wind 'V_LidOpen' ];
end

expRunName = sprintf('%s_Run%d/',expName,runNum);

DataPath = sprintf('%s_Run%d/',[ROOTPath expName '/' expName],runNum);
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_andy/'];
ManualResultsPath = [DataPath 'RESULTS_Manual/'];

runName = sprintf('Run%d',runNum);

load BadPix.mat
CST = load('CST.mat');
load Norm_PIV.mat
%% Test surface-only PIV: Load images and detect surfaces
idx = 480;

PIVAirDir = dir([LoadPath 'PIV Air/' '*.raw']);
PIVWaterDir = dir([LoadPath 'PIV Water/' '*.raw']);
PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV/' '*.raw']);
FI = 0;
LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
image_index = FI+1:2:LI;

pair_index = (image_index(idx)+1)/2;

PIV1Dir_temp = PIVAirDir;
PIV2Dir_temp = PIVWaterDir;
SurfDir_temp = PIVSurf_LFV_Dir;
ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);

ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
ImageNum_Water2 = PIV2Dir_temp(image_index(idx)+1).name(max(strfind(PIV2Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)+1).name)-4);


PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4)

imagename = [LoadPath 'PIVSURF Water/' expName '_' runName '_PIVSURF Water_' ImageNum_Water1 '.raw'];
[IM1] = load_Image_IOCoreView_12MP(imagename);
% Remove Bad Pixels and interpolate
for iiii = 1:length(BadPix_SurfW)
    IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
end
PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;
%-%-%-%-%-%-%-%-%-%-% Load PIVSurf W2
imagename = [LoadPath 'PIVSURF Water/' expName '_' runName '_PIVSURF Water_' ImageNum_Water2 '.raw'];
[IM1] = load_Image_IOCoreView_12MP(imagename);
% Remove Bad Pixels and interpolate
for iiii = 1:length(BadPix_SurfW)
    IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
end
PIVSurf_W2_Raw = fillmissing(IM1,'pchip');
PIVSURF_W2 = PIVSurf_W2_Raw./Norm_PIVSurfW2;

Check_W1 = mean(imbinarize(PIVSURF_W1,30));
Check_W1(1:300) = NaN;
Check_W2 = mean(imbinarize(PIVSURF_W2,30));
Check_W2(1:300) = NaN;

if false && (sum(Check_W1<0.05)>10 || sum(Check_W2<0.05)>10)

    CST.isPIVWater = 0;
    disp(['WINDSHIELD WIPER in WATER Surface!'])
    
    % We skip Water calculations when windshield wipers in the FoV of
    % PIVSurf Water

else
    CST.isPIVWater = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIV images - CONVERT TO MAT
    % First PIV image Water
    filename = [LoadPath 'PIV Water/' expName '_' runName '_PIV Water_' ImageNum_Water1 '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    % Remove Bad Pixels and interpolate
    for iiii = 1:length(BadPix_Water)
        IM1(BadPix_Water(iiii,1),BadPix_Water(iiii,2)) = NaN;
    end
    IM1_W = fillmissing(IM1,'pchip');

    % Second PIV image Water
    filename = [LoadPath 'PIV Water/' expName '_' runName '_PIV Water_' ImageNum_Water2 '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    % Remove Bad Pixels and interpolate
    for iiii = 1:length(BadPix_Water)
        IM1(BadPix_Water(iiii,1),BadPix_Water(iiii,2)) = NaN;
    end
    IM2_W = fillmissing(IM1,'pchip');

    % Camera Angle Correction
    [PIVWater1_CamAngle] = PIVWater_CamAngle_Correction(IM1_W); %Image 1
    [PIVWater2_CamAngle] = PIVWater_CamAngle_Correction(IM2_W); %Image 2
    PIVWater1_CamAngle(PIVWater1_CamAngle<0) = 0;
    PIVWater2_CamAngle(PIVWater2_CamAngle<0) = 0;
    % Correct less brightness in second image
    PIVWater2_CamAngle = PIVWater2_CamAngle*mean(PIVWater1_CamAngle(:),'omitnan')/mean(PIVWater2_CamAngle(:),'omitnan');
    PIVWater2_CamAngle(PIVWater2_CamAngle>1023) = 1023;

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
    %-%-%-%-%-%-%-%-%-%-% Water 1
    % Lens distortion correction
    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % TO BE DONE!!!!! BUT PROBABLY NOT NEEDED

    % Camera Angle Correction
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

    % Extract surface
    [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');

    % Mask PIVWater1
    [Mask1_W] = PIVWater_Mask(PIV1_W, PIVW1_Surface);

    [h, w] = size(PIV1_W);

    %-%-%-%-%-%-%-%-%-%-% Water 2
    % Lens distortion correction
    [PIVSurfW2_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W2);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE
    % DONE!!!!! PROBABLY NOT NEEDED

    % Camera Angle Correction
    [PIVSurfW2_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW2_Undistorted);

    % Extract surface
    [BadFramePIVSurfW2, XPIVSurfW2_Surface, PIVSurfW2_Surface, XPIVW_PIVSurfW2_Surface, PIVW_PIVSurfW2_Surface, PIVW2_Surface,T2_W,RotAngle_W,DY_W] = ExtractSurface_PIVSurfWater(PIVSurfW2_CamAngle,PIV2_W,'2');

    % Mask PIVWater2
    [Mask2_W] = PIVWater_Mask(PIV2_W, PIVW2_Surface);

    %% Compute Surface Velocity
    IntrWndw_surf = [256,128];
    GrdSpc_surf = [128,64];
    
    SurfVel = ComputeSurfaceVelocity(PIV1_W, PIVW1_Surface, PIV2_W, PIVW2_Surface,Mask1_W,Mask2_W,IntrWndw_surf,GrdSpc_surf);
    % PIV1 = PIV1_W;
    % PIV1_Surface = PIVW1_Surface;
    % PIV2 = PIV2_W;
    % PIV2_Surface = PIVW2_Surface;
    % Mask1 = Mask1_W;
    % Mask2 = Mask2_W;
    % IntrWndw = IntrWndw_surf;
    % GrdSpc = GrdSpc_surf;

end

%% Surface PIV

function [SurfVel] = ComputeSurfaceVelocity(PIV1, PIV1_Surface, PIV2, PIV2_Surface, Mask1, Mask2, IntrWndw, GrdSpc)
[h, w] = size(PIV1); % Image height and width
H = 50; %height of near-surface region

[Xp,Zp] = meshgrid(1:w,1:H);

PIV1p = interp2(1:w,(1:h)',PIV1,Xp,Zp+PIV1_Surface-H/2);
PIV2p = interp2(1:w,(1:h)',PIV2,Xp,Zp+PIV2_Surface-H/2);
%% PIV
[CompVelInt] = ComputeVelocities_Surface_Quick_NoFilt_Deform_Water(PIV1p,PIV2p,0*PIV1p + 1,0*PIV2p+1,IntrWndw,GrdSpc);
SurfVel = CompVelInt;
SurfVel.delta_z = SurfVel.delta_z + PIV1_Surface-PIV2_Surface; %double check this
end