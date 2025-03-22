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
        VideoFold = [ROOTPath 'video demonstration\'];
        SaveVideoINTdelx = [VideoFold '\Combined_AirWater\INTdelx\'];
        SaveVideoVorticity = [VideoFold '\Combined_AirWater\Vorticity\'];
        
        % Air
        AirPath = [ResultsPath 'Air\'];
        SavePIVAirPath = [AirPath 'PIV_Velocities_raw\'];
        SaveSurfAirPath = [AirPath 'Surfaces\'];
        FieldsAirPath  = [AirPath 'CALCULATED_FIELDS\'];
        SaveCartAirPath  = [FieldsAirPath 'Cartesian Fields\Velocity\'];
        SaveCartDiffAirPath  = [FieldsAirPath 'Cartesian Fields\Gradients\'];
        SavePressAirPath  = [FieldsAirPath 'Cartesian Fields\Pressure\'];
        SaveTransfoAirPath  = [AirPath 'transfo\'];
        SaveVideoAirINTdelx = [VideoFold '\Air\INTdelx\'];
        SaveVideoAirVort = [VideoFold '\Air\Vorticity\'];
        
        % Water
        WaterPath = [ResultsPath 'Water\'];
        SavePIVWaterPath = [WaterPath 'PIV_Velocities_raw\'];
        SaveSurfWaterPath = [WaterPath 'Surfaces\'];
        FieldsWaterPath  = [WaterPath 'CALCULATED_FIELDS\'];
        SaveCartWaterPath  = [FieldsWaterPath 'Cartesian Fields\Velocity\'];
        SaveCartDiffWaterPath  = [FieldsWaterPath 'Cartesian Fields\Gradients\'];
        SavePressWaterPath  = [FieldsWaterPath 'Cartesian Fields\Pressure\'];
        SaveTransfoWaterPath  = [WaterPath 'transfo\'];
        SaveComboPIVPath = [ResultsPath 'Combined_AirWater\'];
        
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
        
        if ~exist(SaveComboPIVPath, 'dir')
            mkdir(SaveComboPIVPath);
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
        %         save([ResultsPath expName '_' runName '_Parameters.mat'],'CST')
        
        %% Frame to process
        PIVAirDir = dir([LoadPath 'PIV Air\' '*.raw']);
        PIVWaterDir = dir([LoadPath 'PIV Water\' '*.raw']);
        PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV\' '*.raw']);
        FI = 0;
        LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
        
        %%% Pairs already computed
        Dir = dir([SavePressAirPath '*.mat']);
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
        
        disp(runName)
        
        %%% Pairs already computed
        DirA_Vort = dir([SaveCartDiffAirPath '*.mat']);
        DirW_Vort = dir([SaveCartDiffWaterPath '*.mat']);
        DirA_Vel = dir([SavePIVAirPath '*.mat']);
        DirW_Vel = dir([SavePIVWaterPath '*.mat']);
        
        %% Processing frames
        
        Idx = [];
        for i = 1:15
            Idx = [Idx 249+i:15:numel(DirA_Vort)];
        end
        
        for idx = Idx(2:end) %1:numel(DirA_Vort) %---- MAIN LOOP ----%
            
            PairNum = DirA_Vort(idx).name(end-7:end-4);
            
            %% Combo PIV - Air and Water
            %%% Put together PIV Air and PIV Water
            %%% I will conserve this part but will do the calculation in a
            %%% second time. Now I just need PIV_Surface_PIVAir
            
            % Load surfaces, velocities and gradients
            load([SavePIVAirPath DirA_Vel(idx).name]);
            load([SavePIVWaterPath DirW_Vel(idx).name]);
            load([SaveCartDiffAirPath DirA_Vort(idx).name]);
            load([SaveCartDiffWaterPath DirW_Vort(idx).name]);
            
            % Further roto/translation of water PIV image to match flat surface (from
            % ExtractSurface_PIVSurfWater.m)
            RotAngle_W = (674.955-684.981)/(5250.53+477.054);
            DY_W = 4;
            
            % Tform from PIV Water to PIVSurf Water
            T2_W_from_landmarks_PIVWater_to_PIVSurf_Water
            T2_W = T2; clear T2
            
            % Retrieve XPIVSurfW1_Surface from Main1_PIV
            size2_PIVSurfW1_CamAngle = 4179; % size(PIVSurfW1_CamAngle,2) - look at Main1
            XPIVSurfW1_Surface = 501:size2_PIVSurfW1_CamAngle-40;
            
            % Retrieve size(PIV1_W) and size(PIV1_A) from Main1_PIV
            PIV1_W = zeros(3087,4151);
            PIV1_A = zeros(3088,4143);
            
            % Add vorticity to the matched fields
            CompVelAir.Vorticity = cartDiffAir.Vorticity;
            CompVelWater.Vorticity = cartDiffWater.Vorticity;
            
            % Matching PIVSurf Air with CompVelWater in PIVAir coordinates
            Ix1 = find(PixRes_Air.XPIV_LFV_Surface == -46); % first point of PIVWater in PIVAir coordinates (in pixel resolution)
            IxEnd = find(PixRes_Air.XPIV_LFV_Surface == 4199); % last point of PIVWater in PIVAir coordinates (in pixel resolution)
            PIV_Surface_PIVAir = PixRes_Air.PIV_LFV_Surface(Ix1:IxEnd);
            PIV_Surface_PIVAir = PIV_Surface_PIVAir(CST.IW/2:CST.GS:end-CST.GS)/CST.GS; % PIVSurf Air matching PIVWater (PIVAir resolution)
            Ix = -46:4199;
            xPIV_Surface_PIVAir = Ix(CST.IW/2:CST.GS:end-CST.GS)/CST.GS;
            
            %-%-%-%-%-%-%-%-%-%-% Match PIV Air coordinates
            [CompVelWater_PIVAir,CompVelWater_PIVSurfW,Uinv,Vinv] = Match_PIV_coordinates(RotAngle_W,DY_W,T2_W,PixRes_Water1.XPIVW_PIVSurfW1_Surface,PixRes_Water1.PIVW_PIVSurfW1_Surface,XPIVSurfW1_Surface,PixRes_Air.XPIV_LFV_Surface,PixRes_Air.PIV_LFV_Surface,CompVelWater,PIV1_W,PIV1_A);
            
            % Put together velocity in air (in m/s) and velocity in water (cm/s)
            [Combo_PIV.delta_x] = Combine_PIV(CompVelAir.delta_x*CST.DX/CST.DT, CompVelWater_PIVAir.delta_x*CST.DX/CST.DT_W*10, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            [Combo_PIV.delta_z] = Combine_PIV(CompVelAir.delta_z*CST.DX/CST.DT, CompVelWater_PIVAir.delta_z*CST.DX/CST.DT_W*10, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            [Combo_PIV.INTdelx] = Combine_PIV(CompVelAir.INTdelx*CST.DX/CST.DT, CompVelWater_PIVAir.INTdelx*CST.DX/CST.DT_W*10, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            [Combo_PIV.INTdelz] = Combine_PIV(CompVelAir.INTdelz*CST.DX/CST.DT, CompVelWater_PIVAir.INTdelz*CST.DX/CST.DT_W*10, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            [Combo_PIV.dcor] = Combine_PIV(CompVelAir.dcor, CompVelWater_PIVAir.dcor, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            [Combo_PIV.Vorticity] = Combine_PIV(CompVelAir.Vorticity/CST.DT, CompVelWater_PIVAir.Vorticity/CST.DT_W*10, xPIV_Surface_PIVAir, PIV_Surface_PIVAir);
            CST.isCombo = 1;
            Combo_PIV.xPIV_Surface_PIVAir = xPIV_Surface_PIVAir;
            Combo_PIV.PIV_Surface_PIVAir = PIV_Surface_PIVAir;
            
            %% Save Pictures
            FIELDS = {'Vorticity','INTdelx'};
            % FIELD = 'Vorticity'; % NOTE: INTdelx in m/s x 10^{-1) in water, in m/s in air
            % FIELD = 'INTdelx'; % NOTE: INTdelx in cm/s in water, in m/s in air
            for j = 1:2
                FIELD = FIELDS{j};
                %%% Build colormap with white equivalent to the zero of the field
                CMap(1,:) = [ 0 0 0 ]; % black
                CMap(2,:) = [ 0.5 0 0.75 ]; % dark purple
                CMap(3,:) = [ 0 0 1 ]; % blue
                CMap(4,:) = [ 0.3333 1 1 ]; % cyan
                CMap(5,:) = [ 1 1 1 ]; % white
                CMap(6,:) = [ 0 1 0 ]; % green
                CMap(7,:) = [ 1 1 0 ]; % green-yellow
                CMap(8,:) = [ 1 0.5 0 ]; % yellow
                CMap(9,:) = [ 1 0 0]; % red
                Step = 50*size(CMap,1);
                R = interp1(1:size(CMap,1),CMap(:,1)',1:(size(CMap,1)-1)/Step:size(CMap,1),'linear')';
                G = interp1(1:size(CMap,1),CMap(:,2)',1:(size(CMap,1)-1)/Step:size(CMap,1),'linear')';
                B = interp1(1:size(CMap,1),CMap(:,3)',1:(size(CMap,1)-1)/Step:size(CMap,1),'linear')';
                Colormap = [R G B];
                
                figure;
                Ximg = (0:size(Combo_PIV.(FIELD),2)-1)*CST.DX;
                Zimg = fliplr(((0:size(Combo_PIV.(FIELD),1)-1)-(size(Combo_PIV.(FIELD),1)+1-mean(Combo_PIV.PIV_Surface_PIVAir)))*CST.DX);
                imagesc(Ximg,Zimg,Combo_PIV.(FIELD))
                %     eval(['imagesc(Combo_PIV.' FIELD ';'])
                set(gca,'YDir','normal')
                hold on
                plot(Ximg,((size(Combo_PIV.(FIELD),1)-Combo_PIV.PIV_Surface_PIVAir)-mean(size(Combo_PIV.(FIELD),1)+1-Combo_PIV.PIV_Surface_PIVAir))*CST.DX,'Color',[0.5 0.5 0.5],'LineWidth',1.25)
                ElapsedTime = round(str2double(PairNum)/15,2);
                title(['Elapsed time: ' num2str(ElapsedTime) ' s']);
                %     CAXIS = [-abs(max(max(Combo_PIV.Vorticity)))*0.9 abs(max(max(Combo_PIV.Vorticity)))*0.9 ]; % for Vorticity
                %     CAXIS = [-abs(max(max(Combo_PIV.INTdelx)))*0.9 abs(max(max(Combo_PIV.INTdelx)))*0.9 ]; % for INTdelx
                CAXIS = [-abs(max(max(Combo_PIV.(FIELD))))*0.9 abs(max(max(Combo_PIV.(FIELD))))*0.9]; % for general field
                caxis(CAXIS)
                xlabel('x (m)');
                ylabel('z (m)');
                cb = colorbar;
                cb.Label.String = '\omega_a (s^{-1}),    \omega_w (s^{-1} x 10^{-1})'; %for Vorticity
                if j == 2
                    cb.Label.String = 'u_a (ms^{-1}),    u_w (ms^{-1} x 10^{-1})'; % for INTdelx
                end
                cb.Label.Rotation = 270;
                cb.Label.FontSize = 12;
                cb.Label.VerticalAlignment = "bottom";
                cb.FontName = 'Times New Roman';
                colormap(gca,Colormap);
                if j == 1
                    SaveDir = SaveVideoVorticity;
                    print(gcf,[SaveVideoVorticity expName '_' PairNum '_VORTICITY.png'],'-dpng','-r600'); %for Vorticity
                else
                    SaveDir = SaveVideoINTdelx;
                    print(gcf,[SaveVideoINTdelx expName '_' PairNum '_INTDELX.png'],'-dpng','-r600') % for INTdelx
                end
            end
            close all
            
            %% Screen output
            disp(['Pair ' PairNum ' done.']);
            toc
        end
        
    end
    
end