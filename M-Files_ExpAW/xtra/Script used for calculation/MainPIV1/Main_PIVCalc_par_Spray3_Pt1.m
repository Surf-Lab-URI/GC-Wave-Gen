% This MATLAB script loads raw data, extracts surface and computes PIV velocities.
tic
clear
close all
clc

ROOTPath = 'D:\data\EXPERIMENTS\'; % Spray4

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

for i = 3%2%:length(ExpDir) % Loop on the number of experiments

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
    
    for ii = 1%[1,2,9:16]% 1:numel(RunDir) % Loop on the number of Runs

        runName = RunDir(ii).name(20:end);
        if strcmp(ExpAW,'6')
            runName = RunDir(ii).name(28:end);
        end

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

        %%% Save Parameters
        save([ResultsPath expName '_' runName '_Parameters.mat'],'CST')

        %% Frame to process
        PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
        PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
        PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV\' '*.raw']);
        FI = 0;
        LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);

        %%% Pairs already computed
        Dir = dir([SavePressAirPath '*.mat']);
        Dir = [];
        Pair_Done = [];
        try
            if numel(PIVAirDir)<=10000
                for iii = 1:length(Dir)
                    Pair_Done = [Pair_Done {Dir(iii).name(end-7:end-4)}];
                end
            else
                Pair_Done = [Pair_Done {Dir(iii).name(end-8:end-4)}];
            end
        catch
        end
        %% Processing frames
        image_index = FI+1:2:LI;

        for idx = 1000%1:750%numel(image_index) % Main Loop

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
            imagename = [LoadPath 'PIVSurf Air - LFV\' expName '_' runName '_PIVSurf Air - LFV_' PairNum '.raw'];
            [IM1] = (load_Image_IOCoreView_48MP(imagename));
            PIVSurf_A_Raw = IM1;
            PIVSURF_A = PIVSurf_A_Raw./Norm_LFV;

            % If the mean of imbinarized PIVSurf Air - LFV goes to zero at some point, it means
            % there's a drop due to the windshield wiper and we skip that image
            Check = mean(imbinarize(PIVSURF_A,100));
            if sum(Check<0.01)>10

                CST.isPIVAir = 0;
                disp(['WINDSHIELD WIPER in AIR Surface!'])
                continue
                % We skip Air calculations when windshield wipers in the FoV of
                % PIVSurf Air - LFV

            else
                CST.isPIVAir = 1;

                % Load First PIV image Air
                filename = [LoadPath 'PIV Air\' expName '_' runName '_PIV Air_' ImageNum_Air1 '.raw'];
                [IM1] = (load_Image_IOCoreView_12MP(filename));
                % Remove Bad Pixels and interpolate
                for iiii = 1:length(BadPix_Air)
                    IM1(BadPix_Air(iiii,1),BadPix_Air(iiii,2)) = NaN;
                end
                IM1_A = fillmissing(IM1,'pchip');

                % Load Second PIV image Air
                filename = [LoadPath 'PIV Air\' expName '_' runName '_PIV Air_' ImageNum_Air2 '.raw'];
                [IM1] = (load_Image_IOCoreView_12MP(filename));
                % Remove Bad Pixels and interpolate
                for iiii = 1:length(BadPix_Air)
                    IM1(BadPix_Air(iiii,1),BadPix_Air(iiii,2)) = NaN;
                end
                IM2_A = fillmissing(IM1,'pchip');

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

                % Lens distortion correction
                [PIVSurfA_Undistorted] = PIVSurfA_LFV_LensDistCorr(PIVSURF_A,cameraParams);
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

            end

            %% Water

            %%%% Check if the windshield wipers in the FoV
            %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W1
            imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water1 '.raw'];
            [IM1] = load_Image_IOCoreView_12MP(imagename);
            % Remove Bad Pixels and interpolate
            for iiii = 1:length(BadPix_SurfW)
                IM1(BadPix_SurfW(iiii,1),BadPix_SurfW(iiii,2)) = NaN;
            end
            PIVSurf_W1_Raw = fillmissing(IM1,'pchip');
            PIVSURF_W1 = PIVSurf_W1_Raw./Norm_PIVSurfW1;
            %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W2
            imagename = [LoadPath 'PIVSurf Water\' expName '_' runName '_PIVSurf Water_' ImageNum_Water2 '.raw'];
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

            if sum(Check_W1<0.05)>10 || sum(Check_W2<0.05)>10

                CST.isPIVWater = 0;
                disp(['WINDSHIELD WIPER in WATER Surface!'])
                continue
                % We skip Water calculations when windshield wipers in the FoV of
                % PIVSurf Water

            else
                CST.isPIVWater = 1;

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
                filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water2 '.raw'];
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

                %%%%%% COORDINATE TRANSFORMATION DATA
                [h, w] = size(PIV1_W);
                [Uorb_W1,Vorb_W1] = Compute_OrbVel_Water(XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, CST, h, w, PIV1_W);

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

            end


            %% PIV CALCULATIONS

            %%%% PIV CALCULATION IN THE AIR!

            if CST.isPIVAir
                if BadFramePIVSurfLFV == 1

                    disp(['NO PIV AIR. Bad PIVSurfAir - LFV. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface detection in the Air
                    continue
                else

                    [CompVelAir] =  PIV_FAB7_threshold_correlation_and_phase_CD_Deform_Fabio(PIV1_A, PIV2_A, Mask_A, IntrWndw_A, GrdSpc_A);

                end
            end

            %%%% PIV CALCULATION IN THE WATER!
            if CST.isPIVWater
                if BadFramePIVSurfW1 == 1 || BadFramePIVSurfW2 == 1

                    if BadFramePIVSurfW1 == 1
                        disp(['NO PIV WATER. Bad PIVSurfW1. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface in Water1
                    elseif BadFramePIVSurfW2 == 1
                        disp(['NO PIV WATER. Bad PIVSurfW2. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface in Water2
                    end
                    continue

                else

                    [CompVelWater] =  ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIV1_W, PIV2_W, Mask1_W, Mask2_W, Mask1_W, IntrWndw_W, GrdSpc_W, Uorb_W1,Vorb_W1);

                end
            end

            %% Combo PIV - Air and Water
            %%% Put together PIV Air and PIV Water
            %%% I will conserve this part but will do the calculation in a
            %%% second time. Now I just need PIV_Surface_PIVAir
            try
                if CST.isPIVAir || CST.isPIVWater
                    % Matching PIVSurf Air with CompVelWater in PIVAir coordinates
                    Ix1 = find(XPIV_LFV_Surface == -46); % first point of PIVWater in PIVAir coordinates (in pixel resolution)
                    IxEnd = find(XPIV_LFV_Surface == 4199); % last point of PIVWater in PIVAir coordinates (in pixel resolution)
                    PIV_Surface_PIVAir = PIV_LFV_Surface(Ix1:IxEnd);
                    PIV_Surface_PIVAir = PIV_Surface_PIVAir(CST.IW/2:CST.GS:end-CST.GS)/CST.GS; % PIVSurf Air matching PIVWater (PIVAir resolution)
                    Ix = -46:4199;
                    xPIV_Surface_PIVAir = Ix(CST.IW/2:CST.GS:end-CST.GS)/CST.GS;

                    %-%-%-%-%-%-%-%-%-%-% Match PIV Air coordinates
                    % % % [CompVelWater_PIVAir,CompVelWater_PIVSurfW,Uinv,Vinv] = Match_PIV_coordinates(RotAngle_W,DY_W,T2_W,XPIVW_PIVSurfW1_Surface,PIVW_PIVSurfW1_Surface,XPIVSurfW1_Surface,XPIV_LFV_Surface,PIV_LFV_Surface,CompVelWater,PIV1_W,PIV1_A);

                    % % % % Put together velocity in air (in m/s) and velocity in water (cm/s)
                    % % % [Combo_PIV.delta_x] = Combine_PIV(CompVelAir.delta_x*CST.DX/CST.DT, CompVelWater_PIVAir.delta_x*CST.DX/CST.DT_W*100, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
                    % % % [Combo_PIV.delta_z] = Combine_PIV(CompVelAir.delta_z*CST.DX/CST.DT, CompVelWater_PIVAir.delta_z*CST.DX/CST.DT_W*100, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
                    % % % [Combo_PIV.INTdelx] = Combine_PIV(CompVelAir.INTdelx*CST.DX/CST.DT, CompVelWater_PIVAir.INTdelx*CST.DX/CST.DT_W*100, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
                    % % % [Combo_PIV.INTdelz] = Combine_PIV(CompVelAir.INTdelz*CST.DX/CST.DT, CompVelWater_PIVAir.INTdelz*CST.DX/CST.DT_W*100, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
                    % % % [Combo_PIV.dcor] = Combine_PIV(CompVelAir.dcor, CompVelWater_PIVAir.dcor, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
                    CST.isCombo = 1;
                end
            catch
                CST.isCombo = 0;
            end

            %% Saving SURFACES
            %%%% Air
            try
                FiltLength = 1000;

                % "Pixel" Resolution
                Build_PixRes_Air

                % "PIV" Resolution
                Build_PIVRes_Air

                %%% Save surface
                SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];
                save([SaveSurfAirPath SurfFileName '.mat'] , 'PIVRes_Air', 'PixRes_Air','CST');
                CST.isSurfAir = 1;

            catch

                disp(['NO SURFACE AIR AVAILABLE for Pair ' PairNum ' : ']);
                CST.isSurfAir = 0;

            end

            %%%% Water
            try
                %%% Surface 1
                % "Pixel Resolution"
                Build_PixRes_Water1

                % "PIV" Resolution
                Build_PIVRes_Water

                %%% Surface 2
                Build_PixRes_Water2

                save([SaveSurfWaterPath SurfFileName '.mat'] , 'PixRes_Water1', 'PIVRes_Water', 'PixRes_Water2','CST');
                CST.isSurfWater = 1;

            catch

                disp(['NO SURFACE WATER AVAILABLE for Pair ' PairNum ' : ']);
                CST.isSurfWater = 0;

            end

            %% SAVING raw PIV
            %%%% Air
            try
                CompVelAir.Surface = (PixRes_Air.PIV_Surface(PIVRes_Air.xPIV) )/CompVelAir.GS;%Fits perfectly with imagesc(CompVelAir.delta_x)
                CompVelAir.pair_index = pair_index;
                CompVelAir.ImageNum_1 = ImageNum_Air1;
                CompVelAir.ImageNum_2 = ImageNum_Air2;
                CompVelAir.PairNum = PairNum;
                CompVelAir.ExpName = ['ExpAW' ExpAW];
                CompVelAir.Acc = Acc;
                CompVelAir.Wind = Wind;
                CompVelAir.Run = runName;
                CompVelAir.SurfacePhase = PIVRes_Air.Phase;
                PIVAirFileName = [ expName '_' runName '_PIV Air_Velocity_' PairNum];
                save([SavePIVAirPath PIVAirFileName '.mat'], 'CompVelAir','CST');
            catch

            end

            %%%% Water
            try
                CompVelWater.Surface = (PixRes_Water1.PIVW1_Surface(PIVRes_Water.xPIV) )/CompVelWater.GS;%Fits perfectly with imagesc(CompVelWater.delta_x)
                CompVelWater.pair_index = pair_index;
                CompVelWater.ImageNum_1 = ImageNum_Water1;
                CompVelWater.ImageNum_2 = ImageNum_Water2;
                CompVelWater.PairNum = PairNum;
                CompVelWater.ExpName = ['ExpAW' ExpAW];
                CompVelWater.Acc = Acc;
                CompVelWater.Wind = Wind;
                CompVelWater.Run = runName;
                CompVelWater.SurfacePhase = PIVRes_Water.Phase;
                PIVWaterFileName = [ expName '_' runName '_PIV Water_Velocity_' PairNum];
                save([SavePIVWaterPath PIVWaterFileName '.mat'], 'CompVelWater', 'CST');
            catch

            end

            %% COORDINATE TRANSFORMATION

            %%% Air
            try

                [transfo_Air] = GenerateTransfo_Air(PixRes_Air, PIVRes_Air, CST);

                % Save transfo_Air
                TransfoFileName = [expName '_Scene' runName '_transfo_' PairNum];
                save([SaveTransfoAirPath TransfoFileName '.mat'] , 'transfo_Air');
            catch

            end

            %%% Water
            try

                [transfo_Water] = GenerateTransfo_Water(PixRes_Water1, PIVRes_Water, CST);

                % Save transfo_Water
                save([SaveTransfoWaterPath TransfoFileName '.mat'] , 'transfo_Water');
            catch

            end

            %% Smoothing PIV

            %%%% Air
            try
                ThrA = 0.4;
                BadSubpixA = 6;
                [Cartesian_Air] = RemoveOutliers_fabio(CompVelAir,ThrA,BadSubpixA);
                Cartesian_Air.Surface = PIVRes_Air.PIV_Surface;
                Cartesian_Air.pair_index = pair_index;
                Cartesian_Air.ImageNum_1 = ImageNum_Water1;
                Cartesian_Air.ImageNum_2 = ImageNum_Water2;
                Cartesian_Air.PairNum = PairNum;
                Cartesian_Air.ExpName = ['ExpAW' ExpAW];
                Cartesian_Air.Acc = Acc;
                Cartesian_Air.Wind = Wind;
                Cartesian_Air.Run = runName;
                Cartesian_Air.Phase = CompVelAir.SurfacePhase;
                Cartesian_Air.xPIV = CompVelAir.xPIV; % The x coordinates of center of IntrWndws
                Cartesian_Air.zPIV = CompVelAir.zPIV; % The y coordinates of center of IntrWndws
                Cartesian_Air.GS = CompVelAir.GS;
                Cartesian_Air.IW = CompVelAir.IW;
                Cartesian_Air.Slope = transfo_Air.AK_smth(1,:);

                PIVRes_Air.Slope =Cartesian_Air.Slope;
                PIVRes_Air.SlopeH = TransformPhaseVector_decay_hor(Cartesian_Air.Slope, PIVRes_Air, imag(transfo_Air.SU(1,:)));
                PIVRes_Air.PhaseH = TransformPhaseVector_decay_hor(PIVRes_Air.Phase, PIVRes_Air, imag(transfo_Air.SU(1,:)));

                % Save smoothed PIV
                CartFileName = [ expName '_' runName '_CartesianAir_' PairNum];
                save([SaveCartAirPath CartFileName '.mat'], 'Cartesian_Air', 'PixRes_Air', 'PIVRes_Air', 'CST');
            catch

            end
            %%%%

            %%%% Water
            try
                ThrA = 0.2;
                BadSubpixW = 20;
                [Cartesian_Water] = RemoveOutliers_fabio(CompVelWater,ThrA,BadSubpixW);
                Cartesian_Water.Surface = PIVRes_Water.PIVW1_Surface;
                Cartesian_Water.pair_index = pair_index;
                Cartesian_Water.ImageNum_1 = ImageNum_Water1;
                Cartesian_Water.ImageNum_2 = ImageNum_Water2;
                Cartesian_Water.PairNum = PairNum;
                Cartesian_Water.ExpName = ['ExpAW' ExpAW];
                Cartesian_Water.Acc = Acc;
                Cartesian_Water.Wind = Wind;
                Cartesian_Water.Run = runName;
                Cartesian_Water.Phase = CompVelWater.SurfacePhase;
                Cartesian_Water.xPIV = CompVelWater.xPIV; % The x coordinates of center of IntrWndws
                Cartesian_Water.zPIV = CompVelWater.zPIV; % The y coordinates of center of IntrWndws
                Cartesian_Water.GS = CompVelWater.GS;
                Cartesian_Water.IW = CompVelWater.IW;
                Cartesian_Water.Slope = transfo_Water.AK_smth(1,:);

                PIVRes_Water.Slope =Cartesian_Water.Slope;
                PIVRes_Water.SlopeH = TransformPhaseVector_decay_hor(Cartesian_Water.Slope, PIVRes_Water, imag(transfo_Water.SU(1,:)));
                PIVRes_Water.PhaseH = TransformPhaseVector_decay_hor(PIVRes_Water.Phase, PIVRes_Water, imag(transfo_Water.SU(1,:)));

                % Save smoothed PIV
                CartFileName = [ expName '_' runName '_CartesianWater_' PairNum];
                save([SaveCartWaterPath CartFileName '.mat'], 'Cartesian_Water', 'PixRes_Water1', 'PIVRes_Water', 'CST');
            catch

            end
            %%%%%

            %% Cartesian Derivatives

            %%%% Air
            try
                %  du(x,z)/dx, du(x,z)/dz, dw(x,z)/dx, dw(x,z)/dz
                [cartDiffAir.u_x, cartDiffAir.u_z] = csapsDiff(Cartesian_Air.u, 0.005, PIVRes_Air.xPIV, PIVRes_Air.zPIV);
                [cartDiffAir.w_x, cartDiffAir.w_z] = csapsDiff(Cartesian_Air.w, 0.005, PIVRes_Air.xPIV, PIVRes_Air.zPIV);
                cartDiffAir.u_z = -cartDiffAir.u_z;
                cartDiffAir.w_z = -cartDiffAir.w_z;

                %  Vorticity: du(x,z)/dz - dw(x,z)/dx
                cartDiffAir.Vorticity = cartDiffAir.u_z - cartDiffAir.w_x;
                cartDiffAir.Mask = Cartesian_Air.Mask;

                FileName = [ expName '_' runName '_cartDiff_' PairNum];
                save([SaveCartDiffAirPath FileName '.mat'], 'cartDiffAir', 'PixRes_Air', 'PIVRes_Air');
            catch

            end
            %%%%


            %%%% Water
            try
                %  du(x,z)/dx, du(x,z)/dz, dw(x,z)/dx, dw(x,z)/dz
                [cartDiffWater.u_x, cartDiffWater.u_z] = csapsDiff(Cartesian_Water.u, 0.005, PIVRes_Water.xPIV, PIVRes_Water.zPIV);
                [cartDiffWater.w_x, cartDiffWater.w_z] = csapsDiff(Cartesian_Water.w, 0.005, PIVRes_Water.xPIV, PIVRes_Water.zPIV);
                cartDiffWater.u_z = -cartDiffWater.u_z;
                cartDiffWater.w_z = -cartDiffWater.w_z;

                %  Vorticity: du(x,z)/dz - dw(x,z)/dx
                cartDiffWater.Vorticity = cartDiffWater.u_z - cartDiffWater.w_x;
                cartDiffWater.Mask = Cartesian_Water.Mask;

                save([SaveCartDiffWaterPath FileName '.mat'], 'cartDiffAir', 'PixRes_Air', 'PIVRes_Air');
            catch

            end
            %%%%

            %% Pressure Calculation and Laplacian

            %%%% Air
            if CST.isPIVAir
            [LaplA, pressAir] = compute_pressure(CST,PIVRes_Air,transfo_Air,cartDiffAir,GrdSpc_A(end),'AIR');

            FileName = [ expName '_' runName '_pressure_' PairNum];
            save([SavePressAirPath FileName '.mat'], 'pressAir', 'LaplA');
            end
            %%%%

            %%%% STILL WORKING ON: needs to be fixed for BCs at the surface
            %%% (Noye and Morton conditions need to be modified for TOP
            %%% left and right corners in water, instead of BOTTOM that are
            %%% used in the air)
            % % % %%%% Water
            % % % if CST.isPIVWater
            % % % [LaplW, pressWater] = compute_pressure(CST,PIVRes_Water,transfo_Water,cartDiffWater,GrdSpc_W(end),'WATER');
            % % % 
            % % % save([SavePressWaterPath FileName '.mat'], 'pressWater', 'LaplW');
            % % % end
            % % % %%%%
            %%%%%

            %% Screen output
            disp(['Pair ' PairNum ' done.']);
toc
        end

        content = [' EXPERIMENT ' expName ' ' runName ' DONE ' ];
        disp(content)
        psw = fscanf(fopen('C:\Users\Administrator\Desktop\content_Fabio.txt','r'),'%s');
        obj = ['PIV MainCalc velocity AIR-WATER ' runName];

        %     % Email to notify finished experiment
        %     send_email('fabio.addona1@gmail.com',psw,'smtp.gmail.com', obj,content)

        % Insert name script
        Script = 'Fabio1-Main_PIVCalc-par-SPRAY3.m';
        send_email('fabio.addona1@gmail.com',psw,'smtp.gmail.com', obj,['SCRIPT ' Script ' Experiment' expName ' ' runName ' Pt1 concluded!'])

    end

end