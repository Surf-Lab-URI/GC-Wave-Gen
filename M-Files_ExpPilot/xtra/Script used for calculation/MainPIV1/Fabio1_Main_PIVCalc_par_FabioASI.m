% This MATLAB script loads raw data, extracts surface and computes PIV velocities.

clear
close all
clc

ROOTPath = '\\spray3\d\data\ExpPilot\'; % FabioASI

ExpPilots = [ 1 1 2 2 3 3 4 4 4 5 5 ];
Scenes = [ 1 2 1 2 1 2 1 2 3 1 2 ] ;

DeltaT_A = [300 300 300 300 120 120 200 200 200 200 80 ]*1d-6;

for i = 1:length(ExpPilots) % Main Loop
    
    expName = num2str(ExpPilots(i));
    
    sceneName = num2str(Scenes(i));
    
    % Define Path
    DataPath = [ROOTPath 'ExpPilot' expName '\' 'ExpPilot' expName '_Scene' sceneName '\' ];
    LoadPath = [DataPath 'RAW\'];
    RawDataPath = [DataPath 'RAW\'];
    ResultsPath = [DataPath 'RESULTS_fabio_par\'];
    
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
    
    %% Bad Frame, Normalization and lens distortion parameters 
    load cameraParams.mat
    
    %% Parameters    
    IntrWndw_A = [128 64 32 16 8]; %[64 32 16 8];
    GrdSpc_A = [64 32 16 8 4]; %[32 16 8 4];
    IntrWndw_W = [128 64 32 16 8];
    GrdSpc_W = [64 32 16 8 4];
    
    CST.DX = 38.106d-06; % meters per pixel (~38 micron/pix) in PIV Air pixel resolution
    CST.DX_W = 41.782d-6; % meters per pixel (~38 micron/pix) in PIV Water pixel resolution
    CST.DT_A = DeltaT_A(i);
    CST.DT_W = 22.22222d-3;
    CST.GS = GrdSpc_A(end);
    CST.IW = IntrWndw_A(end);
    
    
    %%% Save Parameters
    save([ResultsPath 'ExpPilot' expName '_Scene' sceneName '_Parameters.mat'],'CST')
    
    %% Frame to process
    PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
    PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
    FI = 0;
    LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
    
    %% Processing frames
    image_index = FI+1:2:LI;
    
    for idx = 10:numel(image_index) % Main Loop length(FI:LI)
        
        % Indexes for images
        pair_index = (image_index(idx)+1)/2;
        PIV1Dir_temp = PIVAirDir;
        SurfDir_temp = PIVAirDir;
        ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
        ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
        
        ImageNum_Water1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
        ImageNum_Water2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
        PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4);
        
        if i == 1 || i == 1
            ImageNum_Air1 = PIV1Dir_temp(image_index(idx)-1).name(max(strfind(PIV1Dir_temp(image_index(idx)-1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
            ImageNum_Air2 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
        end
        
        %% Air
        
        % Load First PIV image Air
        filename = [LoadPath 'PIV Air\ExpPilot' expName '_Scene' sceneName '_PIV Air_' ImageNum_Air1 '.raw'];
        [IM1] = (load_Image_IOCoreView_12MP(filename));
        IM1_A = IM1;
        
        % Load Second PIV image Air
        filename = [LoadPath 'PIV Air\ExpPilot' expName '_Scene' sceneName '_PIV Air_' ImageNum_Air2 '.raw'];
        [IM1] = (load_Image_IOCoreView_12MP(filename));
        IM2_A = IM1;
        
        % Camera Angle Correction
        [PIVAir1_CamAngle] = PIVAir_CamAngle_Correction(IM1_A); %Image 1
        [PIVAir2_CamAngle] = PIVAir_CamAngle_Correction(IM2_A); %Image 2
        
        %%% Pre-process
        PIV = PIVAir1_CamAngle;
        [PIV] = Pre_process_PIV_Image_Air(PIV); %pre_proc 1st image
        PIV1_A = PIV;
        PIV = PIVAir2_CamAngle;
        [PIV] = Pre_process_PIV_Image_Air(PIV); %pre_proc 2nd image
        PIV2_A = PIV;
        
        %%% Surface detection
        % Load PIVSurf Air - LFV
        imagename = [LoadPath 'PIVSurf Air - LFV\ExpPilot' expName '_Scene' sceneName' '_PIVSurf Air - LFV_' PairNum '.raw'];
        [IM1] = (load_Image_IOCoreView_48MP(imagename));
        PIVSurf_A_Raw = IM1;
        PIVSURF_A = PIVSurf_A_Raw./(smooth(mean(PIVSurf_A_Raw(2400:end,:)),1000)/max(smooth(mean(PIVSurf_A_Raw(2400:end,:)),1000)))';
        
        % Lens distortion correction
        [PIVSurfA_Undistorted] = PIVSurfA_LFV_LensDistCorr(IM1,cameraParams);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE DONE!!!!!
        
        % Camera Angle Correction
        [PIVSurfA_CamAngle] = PIVSurfAir_LFV_CamAngle_Correction(PIVSurfA_Undistorted);
        
        % Extract surface
        [BadFramePIVSurfLFV,XLFV_Surface,LFV_Surface,XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface] = ExtractSurface_PIVSurfAir_LFV_and_correct_surface(PIVSurfA_CamAngle,PIV1_A,idx);
        
        % Mask PIV
        [Mask_A] = PIVAir_Mask(PIV1_A, PIV_Surface);
%         PIV1_A = PIV1_A(11:3074,41:4066);
%         PIV2_A = PIV2_A(11:3074,41:4066);
%         Mask_A = Mask_A(11:3074,41:4066);
%         %%% NOTE: This last modification after masking must be implemented
%         %%% in all the subsequent steps, and also in the previous (for
%         %%% comparison!!)
%         [XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface] = Modify_after_PIVAir_Mask(XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface);
        
        %%% PIV CALCULATIONS
        %             [CompVelAir] =  ComputeVelocities_Quick_NoFilt_Deform(PIV1_A, PIV2_A, Mask_A, IntrWndw_A, GrdSpc_A);
        [CompVelAir] =  PIV_FAB7_threshold_correlation_and_phase_CD_Deform_Fabio(PIV1_A, PIV2_A, Mask_A, IntrWndw_A, GrdSpc_A);
        
        %% Water
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIV images - CONVERT TO MAT
        % First PIV image Water
        filename = [LoadPath 'PIV Water\ExpPilot' expName '_Scene' sceneName '_PIV Water_' ImageNum_Water1 '.raw'];
        [IM1] = (load_Image_IOCoreView_12MP(filename));
        IM1_W = IM1;
        
        % Second PIV image Water
        filename = [LoadPath 'PIV Water\ExpPilot' expName '_Scene' sceneName '_PIV Water_' ImageNum_Water2 '.raw'];
        [IM1] = (load_Image_IOCoreView_12MP(filename));
        IM2_W = IM1;
        
        % Camera Angle Correction
        [PIVWater1_CamAngle] = PIVWater_CamAngle_Correction(IM1_W); %Image 1
        [PIVWater2_CamAngle] = PIVWater_CamAngle_Correction(IM2_W); %Image 2
        
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
        imagename = [LoadPath 'PIVSurf Water\ExpPilot' expName '_Scene' sceneName' '_PIVSurf Water_' ImageNum_Water1 '.raw'];
        if i == 1
            imagename = [LoadPath 'PIVSurf Water\ExpPilot' expName '_Scene' sceneName' '_PIVSurf Water_' num2str(str2double(ImageNum_Water1)-4,'%.4d') '.raw'];
        end
        [IM1] = load_Image_IOCoreView_12MP(imagename);
        PIVSurf_W1_Raw = IM1;
        PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';
        
        % Lens distortion correction
        [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE DONE!!!!!
        
        % Camera Angle Correction
        [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);
        
        % Extract surface
        [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W);
        
        % Mask PIVWater1
        [Mask1_W] = PIVWater_Mask(PIV1_W, PIVW1_Surface);
        
%-%-%-%-%-%-%-%-%-%-% Load PIVSurf W2
        imagename = [LoadPath 'PIVSurf Water\ExpPilot' expName '_Scene' sceneName' '_PIVSurf Water_' ImageNum_Water2 '.raw'];
        if i == 1
            imagename = [LoadPath 'PIVSurf Water\ExpPilot' expName '_Scene' sceneName' '_PIVSurf Water_' num2str(str2double(ImageNum_Water2)-4,'%.4d') '.raw'];
        end
        [IM1] = load_Image_IOCoreView_12MP(imagename);
        PIVSurf_W2_Raw = IM1;
        PIVSURF_W2 = PIVSurf_W2_Raw./(smooth(mean(PIVSurf_W2_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W2_Raw(2000:end,:)),1000)))';
        
        % Lens distortion correction
        [PIVSurfW2_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W2);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % TO BE DONE!!!!!
        
        % Camera Angle Correction
        [PIVSurfW2_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW2_Undistorted);
        
        % Extract surface
        [BadFramePIVSurfW2, XPIVSurfW2_Surface, PIVSurfW2_Surface, XPIVW_PIVSurfW2_Surface, PIVW_PIVSurfW2_Surface, PIVW2_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW2_CamAngle,PIV2_W);
        
%         % Transform surface to mask PIVWater Surface
%         [XPIVWater2_Surface,PIVWater2_Surface,XPIV_PIVWater2_Surface_CRR,PIV_PIVWater2_Surface_CRR] = Correct_Surface_Water(PIVWater2_CRR,XPIVW_PIVSurfW2_Surface,PIVW_PIVSurfW2_Surface,PIV2_W,idx);
        
        % Mask PIVWater2
        [Mask2_W] = PIVWater_Mask(PIV2_W, PIVW2_Surface);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIV CALCULATIONS
        [CompVelWater] =  ComputeVelocities_Quick_NoFilt_Water(PIV1_W, PIV2_W, Mask1_W, Mask2_W, Mask1_W, IntrWndw_W, GrdSpc_W);
%         [CompVelWater2] =  PIV_FAB7_threshold_correlation_and_phase_CD_Deform_Fabio(PIV1_W, PIV2_W, Mask_W, Mask_W, IntrWndw_W, GrdSpc_W);

%-%-%-%-%-%-%-%-%-%-% Match PIV Air coordinates
        [CompVelWater_PIVAir,CompVelWater_PIVSurfW,Uinv,Vinv] = Match_PIV_coordinates(XPIVW_PIVSurfW1_Surface,PIVW_PIVSurfW1_Surface,XPIVSurfW1_Surface,PIV_LFV_Surface,CompVelWater,PIV1_W,PIV1_A);
        
        %% Put together PIV Air and PIV Water 
        % Matching PIVSurf Air with CompVelWater in PIVAir coordinates
        PIV_Surface_PIVAir = PIV_LFV_Surface(4333:8940);
        PIV_Surface_PIVAir = PIV_Surface_PIVAir(4:4:end)/4;
        
        % Put together velocity in air (in m/s) and velocity in water (cm/s)
        [Combo_PIV.delta_x] = Combine_PIV(CompVelAir.delta_x*CST.DX/CST.DT_A, CompVelWater_PIVAir.delta_x*CST.DX/CST.DT_W*100, PIV_Surface_PIVAir);
        [Combo_PIV.delta_z] = Combine_PIV(CompVelAir.delta_z*CST.DX/CST.DT_A, CompVelWater_PIVAir.delta_z*CST.DX/CST.DT_W*100, PIV_Surface_PIVAir);
        [Combo_PIV.INTdelx] = Combine_PIV(CompVelAir.INTdelx*CST.DX/CST.DT_A, CompVelWater_PIVAir.INTdelx*CST.DX/CST.DT_W*100, PIV_Surface_PIVAir);
        [Combo_PIV.INTdelz] = Combine_PIV(CompVelAir.INTdelz*CST.DX/CST.DT_A, CompVelWater_PIVAir.INTdelz*CST.DX/CST.DT_W*100, PIV_Surface_PIVAir);
        
        %% Saving surfaces
        %%%% Air
        PIVRes_Air.PIV_Surface_PIVAir = PIV_Surface_PIVAir;
        PIVRes_Air.BadFramePIVSurfLFV = BadFramePIVSurfLFV;
        PixRes_Air.BadFramePIVSurfLFV = BadFramePIVSurfLFV;
        PixRes_Air.XLFV_Surface = XLFV_Surface;
        PixRes_Air.LFV_Surface = LFV_Surface;
        PixRes_Air.XPIV_LFV_Surface = XPIV_LFV_Surface;
        PixRes_Air.PIV_LFV_Surface = PIV_LFV_Surface;
        PixRes_Air.PIV_Surface = PIV_Surface;
        
        %%%% Water 
        % Surface 1
        PixRes_Water1.BadFramePIVSurfW1 = BadFramePIVSurfW1;
        PixRes_Water1.XPIVSurfW1_Surface = XPIVSurfW1_Surface;
        PixRes_Water1.PIVSurfW1_Surface = PIVSurfW1_Surface;
        PixRes_Water1.XPIVW_PIVSurfW1_Surface = XPIVW_PIVSurfW1_Surface;
        PixRes_Water1.PIVW_PIVSurfW1_Surface = PIVW_PIVSurfW1_Surface;
        PixRes_Water1.PIVW1_Surface = PIVW1_Surface;
        
        % Surface 2
        PixRes_Water2.BadFramePIVSurfW2 = BadFramePIVSurfW2;
        PixRes_Water2.XPIVSurfW2_Surface = XPIVSurfW2_Surface;
        PixRes_Water2.PIVSurfW2_Surface = PIVSurfW2_Surface;
        PixRes_Water2.XPIVW_PIVSurfW2_Surface = XPIVW_PIVSurfW2_Surface;
        PixRes_Water2.PIVW_PIVSurfW2_Surface = PIVW_PIVSurfW2_Surface;
        PixRes_Water2.PIVW2_Surface = PIVW2_Surface;
        
        %%% Save surface
        SurfFileName = ['ExpPilot' expName '_Scene' sceneName '_Surfaces_' PairNum];
        save([SaveSurfAirPath SurfFileName] , 'PIVRes_Air', 'PixRes_Air','CST');
        save([SaveSurfWaterPath SurfFileName] , 'PixRes_Water1', 'PixRes_Water2','CST');
        
        %% SAVING raw PIV data
        PIVAirFileName = ['ExpPilot' expName '_Scene' sceneName '_PIV Air_Velocity_' PairNum];
        save([SavePIVAirPath PIVAirFileName], 'CompVelAir','Combo_PIV','CST');
        
        PIVWaterFileName = ['ExpPilot' expName '_Scene' sceneName '_PIV Water_Velocity_' PairNum];
        save([SavePIVWaterPath PIVWaterFileName], 'CompVelWater', 'CompVelWater_PIVAir', 'CompVelWater_PIVSurfW','CST');
        
        %% Screen output
        disp(['Pair ' PairNum ' done.']);
        
    end
    
    content = [' EXPERIMENT  ExpPilot' expName ' Scene' sceneName ' DONE ' ];
    disp(content)
    psw = fscanf(fopen('C:\Users\Fabio\Desktop\content_Fabio.txt','r'),'%s');
    obj = 'PIV MainCalc velocity';
    
    %     % Email to notify finished experiment
    %     send_email('fabio.addona1@gmail.com',psw,'smtp.gmail.com', obj,content)
    
    % Insert name script
    Script = 'Fabio1-Main_PIVCalc-par-FabioASI.m';
    send_email('fabio.addona1@gmail.com',psw,'smtp.gmail.com', obj,['SCRIPT ' Script ' ExpPilot' expName ' Scene' sceneName ' concluded!'])
    
    
end