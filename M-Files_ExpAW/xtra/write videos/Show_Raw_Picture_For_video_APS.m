% This MATLAB script loads raw data, extracts surface and computes PIV velocities.
tic
clear
close all
clc

ROOTPath = '\\spray3-10g\d\data\EXPERIMENTS\'; % Spray4

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

for i = 3%:length(ExpDir) % Loop on the number of experiments
    
    ExpAW = ExpDir(i).name(6);
    Acc = ExpDir(i).name(11:14);
    Wind = ExpDir(i).name(17);
    switch Wind
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
    
    RunDir = dir([ExpDir(i).folder '\' ExpDir(i).name '\Exp*Run*']);
    %%% When we want to analyze the "NoFog" runs, we have to modify RunDir.
    %%% But also, we have to modify the script a little bit to skip PIVAir
    %%% analysis
    
    disp(['Exp ' expName ' started!'])
    
    for ii = [1,2,9:16]% 1:numel(RunDir) % Loop on the number of Runs
        
        runName = RunDir(ii).name(20:end);
        if strcmp(ExpAW,'6')
            runName = RunDir(ii).name(28:end);
        end
        
        disp(runName)
        
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
        SaveCartAirPath  = [FieldsAirPath 'Cartesian Fields\Velocity\'];
        SaveCartDiffAirPath  = [FieldsAirPath 'Cartesian Fields\Gradients\'];
        SavePressAirPath  = [FieldsAirPath 'Cartesian Fields\Pressure\'];
        SaveTransfoAirPath  = [AirPath 'transfo\'];
        
        % Water
        WaterPath = [ResultsPath 'Water\'];
        SavePIVWaterPath = [WaterPath 'PIV_Velocities_raw\'];
        SaveSurfWaterPath = [WaterPath 'Surfaces\'];
        FieldsWaterPath  = [WaterPath 'CALCULATED_FIELDS\'];
        SaveCartWaterPath  = [FieldsWaterPath 'Cartesian Fields\Velocity\'];
        SaveCartDiffWaterPath  = [FieldsWaterPath 'Cartesian Fields\Gradients\'];
        SavePressWaterPath  = [FieldsWaterPath 'Cartesian Fields\Pressure\'];
        SaveTransfoWaterPath  = [WaterPath 'transfo\'];
        
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
        
        if ~exist(SaveCartAirPath, 'dir')
            mkdir(SaveCartAirPath);
        end
        if ~exist(SaveCartWaterPath, 'dir')
            mkdir(SaveCartWaterPath);
        end
        
        if ~exist(SaveCartDiffAirPath, 'dir')
            mkdir(SaveCartDiffAirPath);
        end
        
        if ~exist(SaveCartDiffWaterPath, 'dir')
            mkdir(SaveCartDiffWaterPath);
        end
        
        if ~exist(SavePressAirPath, 'dir')
            mkdir(SavePressAirPath);
        end
        
        if ~exist(SavePressWaterPath, 'dir')
            mkdir(SavePressWaterPath);
        end
        
        if ~exist(SaveTransfoAirPath, 'dir')
            mkdir(SaveTransfoAirPath);
        end
        if ~exist(SaveTransfoWaterPath, 'dir')
            mkdir(SaveTransfoWaterPath);
        end
        
        %% Bad Pixels, Normalization and lens distortion parameters
        load cameraParams.mat
        load Norm_PIV.mat
        load BadPix.mat
        
        %% Parameters
        IntrWndw_A = [128 64 32 16 8]; %[64 32 16 8];
        GrdSpc_A = [64 32 16 8 4]; %[32 16 8 4];
        IntrWndw_W = [[256 192 128 96 64 48 32 16 8]*8 32 16 8];%[[256 64 32 16 8]*8 16 8];
        GrdSpc_W = [[128 96 64 48 32 24 16 8 4]*8 16 8 4]; %[[64 32 16 8 4]*8 8 4];
        
        CST.DX = 40.126d-06; % meters per pixel (~40 micron/pix) in PIV Air pixel resolution
        CST.DX_W = 41.242d-6; % meters per pixel (~41 micron/pix) in PIV Water pixel resolution
        CST.DT = DeltaT_A;
        CST.DT_W = 22.22222d-3;
        CST.GS = GrdSpc_A(end);
        CST.IW = IntrWndw_A(end);
        
        % Define physical parameters
        CST.AIR_DENSITY = 1.204;     % air density [kg/m3] at 20°C
        CST.AIR_DVISCOSITY = 1.825e-5;   % dynamic viscosity of air [kg/m*s] at 20°C
        CST.g = 9.81;                 % gravitational acceleration in [m/s2]
        CST.WATER_DEPTH = 0.7;       % Mean water depth in [m]
        CST.WATER_DENSITY = 1000;    % in [kg/m3]
        CST.WATER_DVISCOSITY = 1.0016e-3;   % dynamic viscosity of water [kg/m*s] at 20°C
        CST.SURFACE_TENSION = 0.074; % surface tension in [N/m]
        CST.TOLERANCE = 10e-14;      % numerical tolerance
        
        %% Frame to process
        PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
        PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
        PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV\' '*.raw']);
        PIVSurf_WaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
        FI = 0;
        LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
        
        %%% Pairs already computed
        Dir = [];
        Pair_Done = [];
        
        %% Processing frames
        image_index = FI+1:2:LI;
        
        PN = [];
        PIVA = [];
        PIVW = [];
        Count = 0;
        
        for idx = 250:numel(image_index) % Main Loop
            
            % Indexes for images
            pair_index = (image_index(idx)+1)/2;
            PIV1Dir_temp = PIVAirDir;
            PIV2Dir_temp = PIVWaterDir;
            SurfDir_temp = PIVSurf_LFV_Dir;
            ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
            ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
            
            ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
            ImageNum_Water2 = PIV2Dir_temp(image_index(idx)+1).name(max(strfind(PIV2Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)+1).name)-4);
            
            % For ExpAW1_Run11, PIVSurf Water started 4 images after PIVWater
            if strcmp(ExpAW,'1') && strcmp(runName,'Run11')
                if idx < 5
                    continue
                end
            end
            
            PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4);
            
            %% Continue the paired already analysed
            if ismember(PairNum,Pair_Done)
                continue
            end
            
            disp(['Pair ' PairNum ' STARTED!'])
            
            %% Air
            
            %%%% Check if windshield wipers in the field of view
            % Load PIVSurf Air - LFV
% % %             imagename = [LoadPath 'PIVSurf Air - LFV\' expName '_' runName '_PIVSurf Air - LFV_' PairNum '.raw'];
% % %             [IM1] = (load_Image_IOCoreView_48MP(imagename));
% % %             PIVSurf_A_Raw = IM1;
% % %             PIVSURF_A = PIVSurf_A_Raw./Norm_LFV;
            
            % If the mean of imbinarized PIVSurf Air - LFV goes to zero at some point, it means
            % there's a drop due to the windshield wiper and we skip that image
            
            % Load First PIV image Air
            filename = [LoadPath 'PIV Air\' expName '_' runName '_PIV Air_' ImageNum_Air1 '.raw'];
            [IM1] = (load_Image_IOCoreView_12MP(filename));
            % Remove Bad Pixels and interpolate
            for iiii = 1:length(BadPix_Air)
                IM1(BadPix_Air(iiii,1),BadPix_Air(iiii,2)) = NaN;
            end
            IM1_A = fillmissing(IM1,'pchip');
            
            % Load Second PIV image Air
% % %             filename = [LoadPath 'PIV Air\' expName '_' runName '_PIV Air_' ImageNum_Air2 '.raw'];
% % %             [IM1] = (load_Image_IOCoreView_12MP(filename));
% % %             % Remove Bad Pixels and interpolate
% % %             for iiii = 1:length(BadPix_Air)
% % %                 IM1(BadPix_Air(iiii,1),BadPix_Air(iiii,2)) = NaN;
% % %             end
% % %             IM2_A = fillmissing(IM1,'pchip');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-process
            %%% PIVAir 1
            PIV = IM1_A;
            [PIV] = Pre_process_PIV_Image_Air(PIV); %pre_proc 1st image
            PIV1_A = PIV;
            %%% PIVAir 2
% % %             PIV = IM2_A;
% % %             [PIV] = Pre_process_PIV_Image_Air(PIV); %pre_proc 2nd image
% % %             PIV2_A = PIV;
                        
            %% Water
            
            %%%% Check if the windshield wipers in the FoV
            %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W1
% % %             imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water1 '.raw'];
% % %             [IM1] = load_Image_IOCoreView_12MP(imagename);
% % %             % Remove Bad Pixels and interpolate
% % %             for iiii = 1:length(BadPix_SurfW)
% % %                 IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
% % %             end
% % %             PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
% % %             PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;
            %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W2
% % %             imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water2 '.raw'];
% % %             [IM1] = load_Image_IOCoreView_12MP(imagename);
% % %             % Remove Bad Pixels and interpolate
% % %             for iiii = 1:length(BadPix_SurfW)
% % %                 IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
% % %             end
% % %             PIVSurf_W2_Raw = fillmissing(IM1,'pchip');
% % %             PIVSURF_W2 = PIVSurf_W2_Raw./Norm_PIVSurfW2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIV images - CONVERT TO MAT
            % First PIV image Water
            filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water1 '.raw'];
            [IM1] = (load_Image_IOCoreView_12MP(filename));
            % Remove Bad Pixels and interpolate
            for iiii = 1:length(BadPix_Water)
                IM1(BadPix_Water(iiii,1),BadPix_Water(iiii,2)) = NaN;
            end
            IM1_W = fillmissing(IM1,'pchip');
            
            % Second PIV image Water
% % %             filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water2 '.raw'];
% % %             [IM1] = (load_Image_IOCoreView_12MP(filename));
% % %             % Remove Bad Pixels and interpolate
% % %             for iiii = 1:length(BadPix_Water)
% % %                 IM1(BadPix_Water(iiii,1),BadPix_Water(iiii,2)) = NaN;
% % %             end
% % %             IM2_W = fillmissing(IM1,'pchip');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-process
            %%% PIVWater 1
            PIV = IM1_W;
            [PIV] = Pre_process_PIV_Image_Water_IM1(PIV); %pre_proc 1st image
            PIV1_W = PIV;
            %%% PIVWater 2
% % %             PIV = IM2_W;
% % %             [PIV] = Pre_process_PIV_Image_Water_IM2(PIV); %pre_proc 2nd image
% % %             PIV2_W = PIV;

            %% Gather PIV Air and PIV Water
            PN = [ PN ; PairNum];
            Count = Count+1;
            PIVA(:,:,Count) = PIV1_A;
            PIVW(:,:,Count) = PIV1_W;
            
             %% Screen output
            disp(['Pair ' PairNum ' done.']);
            toc
            
            %%%%%%%% PIV Air
% % %             MeanLevAir = 2145;
% % %             XPIVAir = (1:size(PIV1_A,2))*CST.DX;
% % %             YPIVAir = (size(PIV1_A,1)-(1:size(PIV1_A,1))-(size(PIV1_A,1)-MeanLevAir+1))*CST.DX;
% % %             figure;imagesc(XPIVAir,YPIVAir,PIV2_A);colormap gray;caxis([0 150])
% % %             set(gca,'YDir','normal')
% % %             title(['PIV AIR 1. Elapsed time: ' num2str(round((str2double(PN(i,:))-1)/15,2)) ' s'])
            
            %%%%%%%% PIV Water
% % %             MeanLevWater = 654;
% % %             XPIVWater = (1:size(PIV1_W,2))*CST.DX_W;
% % %             YPIVWater = (size(PIV1_W,1)-(1:size(PIV1_W,1))-(size(PIV1_W,1)-MeanLevWater+1))*CST.DX_W;
% % %             figure;imagesc(XPIVWater,YPIVWater,PIVW(:,:,i));colormap gray;
% % %             set(gca,'YDir','normal')
% % %             title(['PIV WATER 1. Elapsed time: ' num2str(round((str2double(PN(i,:))-1)/15,2)) ' s'])
            
        end
        
    end
    
end