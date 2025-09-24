% This MATLAB script loads raw data, extracts surface and computes PIV velocities.
% tic
clear
% close all
clc

ROOTPath = '/media/surflab/Working24/ExpAW/';

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

i = 5%:length(ExpDir) % Loop on the number of experiments

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

    runName = 'Run2'

    % Define Path
    DataPath = [ROOTPath expName '/' expName '_' runName '/' ];
    LoadPath = [DataPath 'RAW/'];
    RawDataPath = [DataPath 'RAW/'];
    ResultsPath = [DataPath 'RESULTS_andy/'];

    % Air
    AirPath = [ResultsPath 'Air/'];
    SavePIVAirPath = [AirPath 'PIV_Velocities_raw/'];
    SaveSurfAirPath = [AirPath 'Surfaces/'];
    FieldsAirPath  = [AirPath 'CALCULATED_FIELDS/'];
    SaveCartAirPath  = [FieldsAirPath 'Cartesian Fields/Velocity/'];
    SaveCartDiffAirPath  = [FieldsAirPath 'Cartesian Fields/Gradients/'];
    SavePressAirPath  = [FieldsAirPath 'Cartesian Fields/Pressure/'];
    SaveTransfoAirPath  = [AirPath 'transfo/'];

    % Water
    WaterPath = [ResultsPath 'Water/'];
    SavePIVWaterPath = [WaterPath 'PIV_Velocities_raw/'];
    SaveSurfWaterPath = [WaterPath 'Surfaces/'];
    FieldsWaterPath  = [WaterPath 'CALCULATED_FIELDS/'];
    SaveCartWaterPath  = [FieldsWaterPath 'Cartesian Fields/Velocity/'];
    SaveCartDiffWaterPath  = [FieldsWaterPath 'Cartesian Fields/Gradients/'];
    SavePressWaterPath  = [FieldsWaterPath 'Cartesian Fields/Pressure/'];
    SaveTransfoWaterPath  = [WaterPath 'transfo/'];

    %Make directories for results. Uncomment if you want to save the
    %results.
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
    CST.DX = 40.126d-06; % meters per pixel (~40 micron/pix) in PIV Air pixel resolution
    CST.DX_W = 41.242d-6; % meters per pixel (~41 micron/pix) in PIV Water pixel resolution
    CST.DT = DeltaT_A;
    CST.DT_W = 22.22222d-3;
    % CST.GS = GrdSpc_A(end);
    % CST.IW = IntrWndw_A(end);

    %% Frame to process
    PIVAirDir = dir([LoadPath 'PIV Air/' '*.raw']);
    PIVWaterDir = dir([LoadPath 'PIV Water/' '*.raw']);
    PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV/' '*.raw']);
    FI = 0;
    LI = min(length(PIVWaterDir)-1,length(PIVAirDir)-1);
    %% Processing fframes
    image_index = FI+1:2:LI;

    idx = 0433%330%866/2%720/2%numel(image_index) % Main Loop

    % Indexes for images
    pair_index = (image_index(idx)+1)/2;
    PIV1Dir_temp = PIVAirDir;
    PIV2Dir_temp = PIVWaterDir;
    SurfDir_temp = PIVSurf_LFV_Dir;

    ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
    ImageNum_Water2 = PIV2Dir_temp(image_index(idx)+1).name(max(strfind(PIV2Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)+1).name)-4);

    PairNum = SurfDir_temp(pair_index).name(max(strfind(SurfDir_temp(pair_index).name,'_'))+1:length(SurfDir_temp(pair_index).name)-4);
    
    disp(['Pair ' PairNum ' STARTED!'])
    %% Water

    %%%% Check if the windshield wipers in the FoV
    %-%-%-%-%-%-%-%-%-%-% Load PIVSurf W1
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
        disp(['WINDSHIELD WIPER in WATER Surface!'])
    end
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

    %%%%%% COORDINATE TRANSFORMATION DATA
    [h, w] = size(PIV1_W);
    % [Uorb_W1,Vorb_W1] = Compute_OrbVel_Water(XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, CST, h, w, PIV1_W);

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


            % %% PIV CALCULATIONS
            % 
            % %%%% PIV CALCULATION IN THE AIR!
            % 
            % % if CST.isPIVAir
            % %     if BadFramePIVSurfLFV == 1
            % % 
            % %         disp(['NO PIV AIR. Bad PIVSurfAir - LFV. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface detection in the Air
            % %         continue
            % %     else
            % % 
            % %         [CompVelAir] =  PIV_FAB7_threshold_correlation_and_phase_CD_Deform_Fabio(PIV1_A, PIV2_A, Mask_A, IntrWndw_A, GrdSpc_A);
            % % 
            % %     end
            % % end
            % 
            % %%%% PIV CALCULATION IN THE WATER!
            % if CST.isPIVWater
            %     if BadFramePIVSurfW1 == 1 || BadFramePIVSurfW2 == 1
            % 
            %         if BadFramePIVSurfW1 == 1
            %             disp(['NO PIV WATER. Bad PIVSurfW1. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface in Water1
            %         elseif BadFramePIVSurfW2 == 1
            %             disp(['NO PIV WATER. Bad PIVSurfW2. SKIPPED Pair ' PairNum ' : ']); % Reject images with bad surface in Water2
            %         end
            %         continue
            % 
            %     else
            %         % [CompVelWater] =  ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIVWater1_CamAngle, PIVWater1_CamAngle, Mask1_W, Mask2_W, Mask1_W, IntrWndw_W, GrdSpc_W, Uorb_W1,Vorb_W1);
            % 
            %         %[CompVelWater] =  ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIV1_W, PIV2_W, Mask1_W, Mask2_W, Mask1_W, IntrWndw_W, GrdSpc_W, Uorb_W1*0 + 0e-2*CST.DT_W/CST.DX_W,Vorb_W1*0);
            %         [CompVelWater] =  ComputeVelocities_Quick_NoFilt_Deform_Water(PIV1_W, PIV2_W, Mask1_W, Mask2_W, IntrWndw_W, GrdSpc_W);
            % 
            % 
            %         figure(20)
            %         xrange = [3000,3500];
            %         plot(mean(CompVelWater.delta_x(:,xrange(1)/4:xrange(2)/4)*CST.DX_W/CST.DT_W,2,'omitmissing'),(1:size(CompVelWater.delta_x,1))*-4*CST.DX_W,'DisplayName','automated PIV','LineWidth',3)
            %         hold on
            %         %manual results
            %         % x = 2100 to 2600 PairNum 0432 exp5 run 2
            %         % u_man = [1.69e-02, 0.02624490909, 0.02624490909, 0.02624490909, 0.1181020909, 0.1106035455, 0.2005860909, 0.1612187273, 0.02811954545, 0.007498545455];
            %         % z_man = [-3.50E-02,-3.39E-02 -3.41E-02 -3.39E-02 -3.09E-02 -3.10E-02 -2.84E-02 -3.02E-02 -3.33E-02 -3.72E-02];
            % 
            %         % x = 3000 to 3500 PairNum 0259 exp5 run 2
            %         % u_man = [0.0	0	0.01312245455	0.01312245455	0.01312245455	0.009373181818	0.007498545455	0.009373181818	0.01124781818	0.01312245455];
            %         % z_man = [-4.40E-02	-3.57E-02	-3.39E-02	-3.19E-02	-3.10E-02	-3.01E-02	-2.94E-02	-2.88E-02	-2.80E-02	-3.23E-02];
            % 
            %         % x = 3000 to 3500 PairNum 1499 exp5 run 2
            %         % u_man = [9.37E-03	0.01874636364	0.005623909091	0.005623909091	0.007498545455	0.005623909091	0.01124781818	0.009373181818	0.003749272727	0.005623909091];
            %         % z_man = [-1.20E-01	-1.01E-01	-8.95E-02	-8.52E-02	-8.11E-02	-7.70E-02	-6.32E-02	-5.80E-02	-5.50E-02	-5.22E-02];
            % 
            %         % x = 1700 to 2200 PairNum 0432 exp5 run 2
            %         % u_man = [1.31E-02 5.62e-03 5.44E-02 1.78E-01 1.93E-01 2.12E-01 2.01E-01 9.37E-03 2.81E-02];
            %         % z_man = [-3.59E-02 -3.75E-02 -3.21E-02 -2.99E-02 -3.07E-02 -2.86E-02 -3.08E-02 -3.47E-02 -3.36E-02];
            % 
            %         % x = 3000 to 3500 PairNum 0399 exp5 run 2
            %         u_man = [7.50E-03 0.01312245455 0.02437027273 0.1312245455 0.1330991818 0.1349738182 0.1199767273 0.1330991818 0.1349738182 0.04874054545 0.1330991818];
            %         z_man = [-0.035736193 -0.03464328 -0.032766769 -0.030147902 -0.029405546 -0.028539464 -0.029591135 -0.028786916 -0.029343683 -0.030271628 -0.029343683];
            %         plot(u_man,z_man,'.r','DisplayName','manual','MarkerSize',25)
            %         legend('Location', 'southeast','Interpreter','latex')
            %         s = sprintf('%s %s PairNum %s, x = %d to %d pixels (%.2f to %.2f cm)',expName(1:6),runName,PairNum,xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
            %         title(s,'Interpreter','latex')
            %         set(gca,'FontSize',20,'TickLabelInterpreter','latex');
            %         xlabel('u (m/s)','Interpreter','latex');
            %         ylabel('z (m)','Interpreter','latex');
            %         ylim([-0.05,-0.027])
            %         % xlim([0,0.22])
            % 
            %     end
            % end
%% Generate .tif files to feed ML PIV
figure
imagesc(PIV1_W)
colormap gray
hold on
axis equal
ylim([500,800])
%%
ysurf = 665; %677 for 0399, 665 for 0432, %665 for 0259, %800 for 1499
shiftmps = 0;
shift = uint16(round(shiftmps*CST.DT_W/CST.DX_W));

imwrite(uint8(PIV2_W(ysurf:end,:)*0.64.*Mask2_W(ysurf:end,:)),[ImageNum_Water2, '.tif'])
imwrite(uint8(PIV1_W(ysurf:end,:)*0.84.*Mask1_W(ysurf:end,:)),[ImageNum_Water1, '.tif'])
%%
outfname = ['/home/surflab/GitRepos/piv_liteflownet-pytorch/images/demo/DemoOutput/PIV-LiteFlowNet-en/-0_2/flow/', ImageNum_Water1, '_out.flo']; %For Pairnum 0399
[u, v] = read_flo_file(outfname);
figure(20)
hold on
s = sprintf('CNN (piv-liteflownet-pytorch) with %.3f m/s shift',shiftmps);
plot(mean((u+double(shift))*CST.DX_W/CST.DT_W,2,'omitmissing'),(ysurf:(ysurf+size(u,1)-1))'*-1*CST.DX_W,'DisplayName',s,'LineWidth',3)
%% Plot manual and PIV velocity profiles on top of each other for either side of a trough
% figure
xrange = xl;
sauto = sprintf('automated PIV, right of trough (x = %d to %d pixels (%.2f to %.2f cm)',xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
plot(mean(CompVelWater.delta_x(:,xrange(1)/4:xrange(2)/4)*CST.DX_W/CST.DT_W,2,'omitmissing'),(1:size(CompVelWater.delta_x,1))*-4*CST.DX_W,'-r','DisplayName',sauto,'LineWidth',3)
hold on
%manual results
% x = 2100 to 2600 PairNum 0432 exp5 run 2
% u_man = [1.69e-02, 0.02624490909, 0.02624490909, 0.02624490909, 0.1181020909, 0.1106035455, 0.2005860909, 0.1612187273, 0.02811954545, 0.007498545455];
% z_man = [-3.50E-02,-3.39E-02 -3.41E-02 -3.39E-02 -3.09E-02 -3.10E-02 -2.84E-02 -3.02E-02 -3.33E-02 -3.72E-02];

% x = 1700 to 2200 PairNum 0432 exp5 run 2
% u_man = [1.31E-02 5.62e-03 5.44E-02 1.78E-01 1.93E-01 2.12E-01 2.01E-01 9.37E-03 2.81E-02];
% z_man = [-3.59E-02 -3.75E-02 -3.21E-02 -2.99E-02 -3.07E-02 -2.86E-02 -3.08E-02 -3.47E-02 -3.36E-02];


% x = 3000 to 3500 PairNum 0399 exp5 run 2
u_man = [7.50E-03 0.01312245455 0.02437027273 0.1312245455 0.1330991818 0.1349738182 0.1199767273 0.1330991818 0.1349738182 0.04874054545 0.1330991818];
z_man = [-0.035736193 -0.03464328 -0.032766769 -0.030147902 -0.029405546 -0.028539464 -0.029591135 -0.028786916 -0.029343683 -0.030271628 -0.029343683];

sman = sprintf('manual, right of trough (x = %d to %d pixels (%.2f to %.2f cm))',xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
plot(u_man,z_man,'.r','DisplayName',sman,'MarkerSize',25)
legend('Location', 'southeast','Interpreter','latex')
s = sprintf('%s %s PairNum %s',expName(1:6),runName,PairNum);
title(s,'Interpreter','latex')
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
xlabel('u (m/s)','Interpreter','latex');
ylabel('z (m)','Interpreter','latex');
ylim([-0.05,-0.0265])
xlim([0,0.22])
