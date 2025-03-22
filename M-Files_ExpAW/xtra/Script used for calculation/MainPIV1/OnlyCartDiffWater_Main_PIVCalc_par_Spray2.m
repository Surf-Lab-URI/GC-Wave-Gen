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
        
        %% Frame to process
        PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
        PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
        PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV\' '*.raw']);
        FI = 0;
        LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);

        %%% Pairs already computed
        DirW = dir([SaveCartWaterPath '*.mat']);

        disp(runName)
        
        for idx = 1:numel(DirW)% 101:150:701 %1:750%numel(image_index) % Main Loop

            PairNum = DirW(idx).name(end-7:end-4);
            
            %% Cartesian Derivatives

            %%% Load Cartesian in water
            CartFileName = [ expName '_' runName '_CartesianWater_' PairNum];
            load([SaveCartWaterPath CartFileName '.mat']);


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

                FileName = [ expName '_' runName '_cartDiff_' PairNum];
                save([SaveCartDiffWaterPath FileName '.mat'], 'cartDiffWater', 'PixRes_Water1', 'PIVRes_Water');
            catch

            end
            %%%%

            %% Screen output
            disp(['Pair ' PairNum ' done.']);
toc
        end

    end

end