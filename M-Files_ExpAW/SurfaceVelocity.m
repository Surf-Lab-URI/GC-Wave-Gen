clear
clc

ROOTPath = '/media/surflab/Working24/ExpAW/';

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

for i = [1 2 4 5 6 7] % ExpNumber
runNum = 5;

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
load Norm_PIV.mat
%% Test surface-only PIV: Load images and detect surfaces

PIVAirDir = dir([LoadPath 'PIV Air/' '*.raw']);
PIVWaterDir = dir([LoadPath 'PIV Water/' '*.raw']);
% PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV/' '*.raw']);
FI = 0;
LI = length(PIVWaterDir)-1;
image_index = FI+1:2:LI;

parfor idx = 1:length(image_index) 
    try
    CST = load('CST.mat');

    pair_index = (image_index(idx)+1)/2;
    
    % PIV1Dir_temp = PIVAirDir;
    PIV2Dir_temp = PIVWaterDir;
    SurfDir_temp = PIVWaterDir;
    % ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
    % ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
    
    ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
    ImageNum_Water2 = PIV2Dir_temp(image_index(idx)+1).name(max(strfind(PIV2Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)+1).name)-4);
    
    
    PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4)
    % try
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
    
    if (sum(Check_W1<0.05)>10 || sum(Check_W2<0.05)>10)
    
        CST.isPIVWater = 0;
        CST.isSurfWater = 0;
        disp(['WINDSHIELD WIPER in WATER Surface!'])
        
        % Keep track of when the wiper is in the water surface frame
    
    else
        CST.isPIVWater = 1;
        CST.isSurfWater = 1;
    end
    
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
    
    % Camera Angle Correction
    [PIVSurfW2_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW2_Undistorted);
    
    % Extract surface
    [BadFramePIVSurfW2, XPIVSurfW2_Surface, PIVSurfW2_Surface, XPIVW_PIVSurfW2_Surface, PIVW_PIVSurfW2_Surface, PIVW2_Surface,T2_W,RotAngle_W,DY_W] = ExtractSurface_PIVSurfWater(PIVSurfW2_CamAngle,PIV2_W,'2');
    
    % Mask PIVWater2
    [Mask2_W] = PIVWater_Mask(PIV2_W, PIVW2_Surface);
    
    %% Compute Surface Velocity
    IntrWndw_surf = [256,128];
    GrdSpc_surf = [128,64];
    
    H = 50; %Height of near-surface region to use for surface velocity calculation
    SurfVel = ComputeSurfaceVelocity(PIV1_W, PIVW1_Surface, PIV2_W, PIVW2_Surface,Mask1_W,Mask2_W,IntrWndw_surf,GrdSpc_surf,H);
    % PIV1 = PIV1_W;
    % PIV1_Surface = PIVW1_Surface;
    % PIV2 = PIV2_W;
    % PIV2_Surface = PIVW2_Surface;
    % Mask1 = Mask1_W;
    % Mask2 = Mask2_W;
    % IntrWndw = IntrWndw_surf;
    % GrdSpc = GrdSpc_surf;
    
    %% Save Results
    WaterPath = [ResultsPath 'Water/'];
    SaveSurfWaterPath = [WaterPath 'Surfaces/'];
    
    if ~exist(SaveSurfWaterPath, 'dir')
        system(['chmod 777 ' DataPath]);
        mkdir(SaveSurfWaterPath);
    end
    %%%% Water
    PixRes_Water1 = struct();
    % PIVRes_Water = struct();
    PixRes_Water2 = struct();
    
    SurfRes_Water1 = struct();
    SurfRes_Water2 = struct();
    
    
        %%% Surface 1
        % "Pixel Resolution"
        FiltLength = 1000;
        PixRes_Water1.BadFramePIVSurfW1 = BadFramePIVSurfW1;
        PixRes_Water1.XPIVW_PIVSurfW1_Surface = XPIVW_PIVSurfW1_Surface; % PIVWater x-axis in PIV Water coordinates
        PixRes_Water1.PIVW_PIVSurfW1_Surface = PIVW_PIVSurfW1_Surface; % PIVWater surface in PIV Water coordinates
        PixRes_Water1.PIVW_PIVSurfW1_Surface_smth = filtfilt(ones(1,round(FiltLength/3.333))/(round(FiltLength/3.333)), 1, PIVW_PIVSurfW1_Surface); % smoothed version to calculate gradients
        PixRes_Water1.PIVW1_Surface = PIVW1_Surface; % Fits perfectly with imagesc(PIV1_W)
        
        % [XPIVW_LFV_Surface,PIVW_LFV_Surface_smth] = transform_phase_from_PIVAir_to_PIVWater(XPIV_LFV_Surface,PixRes_Air.PIV_LFV_Surface_smth); %%% Find phase from LFV in PIVWater coordinates
        % This is needed to calculate the phase in water space
        % coordinates with long components from LFV
        
        % PixRes_Water1.XLFV_Water = XPIVW_LFV_Surface;
        % PixRes_Water1.LFV_Water_smth = PIVW_LFV_Surface_smth;
        % PixRes_Water1.LFV_Water_smth_phase = angle(hilbert(-PixRes_Water1.LFV_Water_smth+mean(PixRes_Water1.LFV_Water_smth,'omitnan')));
        
        PixRes_Water1.pair_index = pair_index;
        PixRes_Water1.ImageNum_1 = ImageNum_Water1;
        PixRes_Water1.PairNum = PairNum;
        PixRes_Water1.ExpName = ['ExpAW' ExpAW];
        PixRes_Water1.Acc = Acc;
        PixRes_Water1.Wind = Wind;
        PixRes_Water1.Run = runName;
        PixRes_Water1.expRunName = expRunName;
    
        % "PIV" Resolution
        % PIVRes_Water.BadFramePIVSurfW1 = BadFramePIVSurfW1;
        % PIVRes_Water.xPIV = SurfVel.xPIV; % The x coordinates of center of IntrWndws
        % PIVRes_Water.zPIV = SurfVel.zPIV; % The y coordinates of center of IntrWndws
        % % PIVRes_Water.GS = CompVelWater.GS; % Final grid spacing
        % PIVRes_Water.PIVW1_Surface = (PixRes_Water1.PIVW1_Surface(PIVRes_Water.xPIV) )/CompVelWater.GS;%Fits perfectly with imagesc(CompVel.delta_x)
        % PIVRes_Water.pair_index = pair_index;
        % PIVRes_Water.ImageNum_1 = ImageNum_Water1;
        % PIVRes_Water.ImageNum_2 = ImageNum_Water2;
        % PIVRes_Water.PairNum = PairNum;
        % PIVRes_Water.ExpName = ['ExpAW' ExpAW];
        % PIVRes_Water.Acc = Acc;
        % PIVRes_Water.Wind = Wind;
        % PIVRes_Water.Run = runName;
        % PIVRes_Water.PF_Surface = length(PIVRes_Water.zPIV)-PIVRes_Water.PIVW1_Surface+1; % It is needed for transformations;
        % %it's the surface that would be detected on an upside down PIV image . %Fits perfectly with imagesc(flipud(CompVelWater.delta_x))
        % I = ismember(PixRes_Water1.XLFV_Water,PIVRes_Water.xPIV);
        % PIVRes_Water.Phase = PixRes_Water1.LFV_Water_smth_phase(I);
    
        %%% Surface 2
        PixRes_Water2.BadFramePIVSurfW2 = BadFramePIVSurfW2;
        PixRes_Water2.XPIVW_PIVSurfW2_Surface = XPIVW_PIVSurfW2_Surface;
        PixRes_Water2.PIVW_PIVSurfW2_Surface = PIVW_PIVSurfW2_Surface;
        PixRes_Water2.PIVW2_Surface = PIVW2_Surface;
        
        PixRes_Water2.pair_index = pair_index;
        PixRes_Water2.ImageNum_2 = ImageNum_Water2;
        PixRes_Water2.PairNum = PairNum;
        PixRes_Water2.ExpName = ['ExpAW' ExpAW];
        PixRes_Water2.Acc = Acc;
        PixRes_Water2.Wind = Wind;
        PixRes_Water2.Run = runName;
        PixRes_Water2.expRunName = expRunName;
    
        % store surface image resolution surfaces in SurfRes_Water1
        SurfRes_Water1.BadFramePIVSurfW1 = BadFramePIVSurfW1;
        SurfRes_Water1.XPIVSurfW1_Surface = XPIVSurfW1_Surface;
        SurfRes_Water1.PIVSurfW1_Surface = PIVSurfW1_Surface;
        SurfRes_Water1.pair_index = pair_index;
        SurfRes_Water1.ImageNum_1 = ImageNum_Water1;
        SurfRes_Water1.PairNum = PairNum;
        SurfRes_Water1.ExpName = ['ExpAW' ExpAW];
        SurfRes_Water1.Acc = Acc;
        SurfRes_Water1.Wind = Wind;
        SurfRes_Water1.Run = runName;
        SurfRes_Water1.expRunName = expRunName;
    
        % store surface image resolution surfaces in SurfRes_Water2
        SurfRes_Water2.BadFramePIVSurfW2 = BadFramePIVSurfW2;
        SurfRes_Water2.XPIVSurfW2_Surface = XPIVSurfW2_Surface;
        SurfRes_Water2.PIVSurfW2_Surface = PIVSurfW2_Surface;
        SurfRes_Water2.pair_index = pair_index;
        SurfRes_Water2.ImageNum_2 = ImageNum_Water2;
        SurfRes_Water2.PairNum = PairNum;
        SurfRes_Water2.ExpName = ['ExpAW' ExpAW];
        SurfRes_Water2.Acc = Acc;
        SurfRes_Water2.Wind = Wind;
        SurfRes_Water2.Run = runName;
        SurfRes_Water2.expRunName = expRunName;
    
        SurfVel.pair_index = pair_index;
        SurfVel.PairNum = PairNum;
        SurfVel.ExpName = ['ExpAW' ExpAW];
        SurfVel.Acc = Acc;
        SurfVel.Wind = Wind;
        SurfVel.Run = runName;
        SurfVel.expRunName = expRunName;
    
        SavedSurfsWater = struct();
        SavedSurfsWater.CST = CST;
        SavedSurfsWater.PixRes_Water1 = PixRes_Water1;
        SavedSurfsWater.PixRes_Water2 = PixRes_Water2;
        SavedSurfsWater.SurfRes_Water1 = SurfRes_Water1;
        SavedSurfsWater.SurfRes_Water2 = SurfRes_Water2;
        % SavedSurfsWater.PIVRes_Water1 = PIVRes_Water1
        SavedSurfsWater.SurfVel = SurfVel;
    
        SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];
        save([SaveSurfWaterPath SurfFileName '.mat'] , '-fromstruct',SavedSurfsWater);
        disp(['SAVED SURFACE for Pair ' PairNum]);
    
    
    catch

        disp(['NO SURFACE WATER AVAILABLE for Pair ' PairNum ' : ']);
        CST.isSurfWater = 0;

    end

end
end

%% Surface PIV

function [SurfVel] = ComputeSurfaceVelocity(PIV1, PIV1_Surface, PIV2, PIV2_Surface, Mask1, Mask2, IntrWndw, GrdSpc, H)
%% Warp the image so that the surface is a straight line.
%H = 5 recommended
[h, w] = size(PIV1); % Image height and width

[Xp,Zp] = meshgrid(1:w,1:H);

PIV1p = interp2(1:w,(1:h)',PIV1,Xp,Zp+PIV1_Surface-H/2);
PIV2p = interp2(1:w,(1:h)',PIV2,Xp,Zp+PIV2_Surface-H/2);
%% PIV
[CompVelInt] = ComputeVelocities_Surface_Quick_NoFilt_Deform_Water(PIV1p,PIV2p,0*PIV1p + 1,0*PIV2p+1,IntrWndw,GrdSpc);
SurfVel = CompVelInt;
for j = 1:size(SurfVel.delta_z,1)
    SurfVel.delta_z(j,:) = SurfVel.delta_z(j,:) + interp1(1:w,PIV1_Surface-PIV2_Surface,SurfVel.xPIV); %double check this
end
SurfVel.delta_z1 = SurfVel.delta_z1 + PIV1_Surface-PIV2_Surface;

SurfVel.H = H;
SurfVel.PIV1_Surface = PIV1_Surface;
SurfVel.PIV2_Surface = PIV2_Surface;
SurfVel.IntrWndw = IntrWndw;
SurfVel.GrdSpc = GrdSpc;
end