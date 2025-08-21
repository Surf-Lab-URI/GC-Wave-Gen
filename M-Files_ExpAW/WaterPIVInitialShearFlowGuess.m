%% Define Paths
clear
clc
% close all

ROOTPath = '/media/surflab/Working24/ExpAW/';

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

for i = [1 2 3 4 5 6 7] % ExpNumber
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
WaterPath = [ResultsPath 'Water/'];
SaveSurfWaterPath = [WaterPath 'Surfaces/'];

runName = sprintf('Run%d',runNum);

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']);

load BadPix.mat
load Norm_PIV.mat
%% Parameters

IntrWndw_W = [512 256 128 64 32];
GrdSpc_W = [256 128 64 32 16];
% 
% IntrWndw_W_i = [256 128 64 32 16 8];
% GrdSpc_W_i = [128 64 32 16 8 4];

IntrWndw_W_i = [512 256, 128];
GrdSpc_W_i = [256 128, 64];

CST = struct();

CST.DX = 40.126d-06; % meters per pixel (~40 micron/pix) in PIV Air pixel resolution
CST.DX_W = 41.242d-6; % meters per pixel (~41 micron/pix) in PIV Water pixel resolution
CST.DT = DeltaT_A;
CST.DT_W = 22.22222d-3;
% CST.GS = GrdSpc_A(end);
% CST.IW = IntrWndw_A(end);
CST.GS_W = GrdSpc_W;
CST.IW_W = IntrWndw_W;

% Define physical parameters
CST.AIR_DENSITY = 1.204;     % air density [kg/m3] at 20°C
CST.AIR_DVISCOSITY = 1.825e-5;   % dynamic viscosity of air [kg/m*s] at 20°C
CST.g = 9.81;                 % gravitational acceleration in [m/s2]
CST.WATER_DEPTH = 0.7;       % Mean water depth in [m]
CST.WATER_DENSITY = 1000;    % in [kg/m3]
CST.WATER_DVISCOSITY = 1.0016e-3;   % dynamic viscosity of water [kg/m*s] at 20°C
CST.SURFACE_TENSION = 0.074; % surface tension in [N/m]
CST.TOLERANCE = 10e-14;      % numerical tolerance
CST.isPIVWater = 0;
CST.isSurfWater = 0;

%Uncomment to save parameters used for processing
%%% Save Parameters
save([ResultsPath expName '_' runName '_Parameters.mat'],'CST')
%% Loop through pairs
PairNumsInt = 0:floor((length(PIVWaterDir)+1)/2)-1;
numPairs = length(PairNumsInt);

% PIVAirDir = dir([LoadPath 'PIV Air/' '*.raw']);
PIVWaterDir = dir([LoadPath 'PIV Water/' '*.raw']);
% PIVSurf_LFV_Dir = dir([LoadPath 'PIVSurf Air - LFV/' '*.raw']);
FI = 0;
LI = length(PIVWaterDir)-1;
image_index = FI+1:2:LI;


 parfor idx = 1:length(image_index)
    PairNumInt = PairNumsInt(idx);
    PairNum = sprintf('%04d',PairNumInt)

    %% Load surface info
    SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];
    PixRes_Water1 = struct();
    PixRes_Water2 = struct();
    SurfRes_Water1 = struct();
    SurfRes_Water2 = struct();
    SurfVel = struct();
    % CST = struct();
    try
        SavedSurfsWater = load([SaveSurfWaterPath SurfFileName '.mat']);
        PixRes_Water1 = SavedSurfsWater.PixRes_Water1;
        PixRes_Water2 = SavedSurfsWater.PixRes_Water2;
        SurfRes_Water1 = SavedSurfsWater.SurfRes_Water1;
        SurfRes_Water2 = SavedSurfsWater.SurfRes_Water2;
        SurfVel = SavedSurfsWater.SurfVel;
        % CST = SavedSurfsWater.CST;
        % CST.isPIVWater = SavedSurfsWater.isPIVWater;
        isSurfWater = SavedSurfsWater.CST.isSurfWater;
        isPIVWater = SavedSurfsWater.CST.isPIVWater;
        surfVelMean = mean(SurfVel.delta_x,'all','omitmissing');%*SavedSurfsWater.CST.DX_W/SavedSurfsWater.CST.DT_W;
        XPIVW_PIVSurfW1_Surface = PixRes_Water1.XPIVW_PIVSurfW1_Surface;
        PIVW_PIVSurfW1_Surface = PixRes_Water1.PIVW_PIVSurfW1_Surface;
        PIVW1_Surface= PixRes_Water1.PIVW1_Surface;
        XPIVW_PIVSurfW2_Surface = PixRes_Water2.XPIVW_PIVSurfW2_Surface;
        PIVW_PIVSurfW2_Surface = PixRes_Water2.PIVW_PIVSurfW2_Surface;
        PIVW2_Surface= PixRes_Water2.PIVW2_Surface;
        disp(['Loaded surfaces from ' SaveSurfWaterPath SurfFileName '.mat'])
    
        %% Load and Preprocess PIV Water images
        pair_index = (image_index(idx)+1)/2;
        
        PIV2Dir_temp = PIVWaterDir;
    
        ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
        ImageNum_Water2 = PIV2Dir_temp(image_index(idx)+1).name(max(strfind(PIV2Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)+1).name)-4);
    
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
                
        % Mask PIVWater1
        [Mask1_W] = PIVWater_Mask(PIV1_W, PIVW1_Surface-20);
        
        [h, w] = size(PIV1_W);
                
        % Mask PIVWater2
        [Mask2_W] = PIVWater_Mask(PIV2_W, PIVW2_Surface-20);
        % CompVelWater = ComputeVelocities_Quick_NoFilt_Deform_Water(PIV1_W, PIV2_W, Mask1_W, Mask2_W, IntrWndw_W_i, GrdSpc_W_i);
        CompVelWater_i = ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIV1_W,PIV2_W,Mask1_W,Mask2_W,Mask1_W,IntrWndw_W_i,GrdSpc_W_i,0*Mask1_W,0*Mask1_W,true);
        deepVelMean = mean(CompVelWater_i.delta_x(floor(size(CompVelWater_i.delta_x,1)*0.5):end,:),'all','omitmissing');
        
        u_bar_i = mean(CompVelWater_i.delta_x,2,'omitmissing');

        shearMask = ((u_bar_i - deepVelMean) > 0.05*(surfVelMean-deepVelMean));
        shearMaskIs = find(shearMask == 1);
        bottomShearLayerI = shearMaskIs(end);
        bottomShearLayer = CompVelWater_i.zPIV(bottomShearLayerI);
        shearGuess = zeros(h,w);
        mwl = floor(mean(PIVW1_Surface,'all','omitmissing'));
        dudy = (deepVelMean-surfVelMean)/(bottomShearLayer-mwl);
        [X1,Y1] = meshgrid(1:w, 1:h);
        shearGuess(bottomShearLayer:end,:) = deepVelMean;
        shearGuess(1:bottomShearLayer,:) = (Y1(1:bottomShearLayer,:)-bottomShearLayer)*dudy + deepVelMean;
        shearGuess(1:mwl,:) = surfVelMean;
        shearGuess = movmean(shearGuess,50,1,'omitmissing');
        SavedSurfsWater.SurfVel.mwl = mwl;
        SavedSurfsWater.SurfVel.bottomShearLayer = bottomShearLayer;
        SavedSurfsWater.SurfVel.dudy = dudy;

        save([SaveSurfWaterPath SurfFileName '.mat'] , '-fromstruct',SavedSurfsWater);
        disp(['SAVED SURFACE for Pair ' PairNum]);
        % th = 170;
        % PIV1 = PIV1_W;
        % PIV2 = PIV2_W;
        % 
        % PIV1(PIV1 < th) = 0;
        % PIV2(PIV2*.9 < th) = 0;
        % CompVelWater = ComputeVelocities_Quick_NoFilt_Deform_Water_InitVel(PIV1,PIV2,Mask1_W,Mask2_W,IntrWndw_W,GrdSpc_W,shearGuess,zeros(h,w),false);
        % CompVelWater = ComputeVelocities_Quick_NoFilt_Deform_Water_InitVel(PIVWater1_CamAngle,PIVWater2_CamAngle,Mask1_W,Mask2_W,IntrWndw_W, GrdSpc_W,zeros(h,w),zeros(h,w),false); 

        % CompVelWater = ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIV1,PIV2,Mask1_W,Mask2_W,Mask1_W,IntrWndw_W, GrdSpc_W,shearGuess,zeros(h,w),false);
        % CompVelWater = ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(PIVWater1_CamAngle,PIVWater2_CamAngle,Mask1_W,Mask2_W,Mask1_W,IntrWndw_W, GrdSpc_W,zeros(h,w),zeros(h,w),false); 
    catch
        SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];
        disp(['Failed on file ' SaveSurfWaterPath SurfFileName '.mat'])
    end


 end  
 end
%%    
        figure(20)
            xrange = [3000,3500];
            plot(mean(CompVelWater_i.delta_x(:,floor(xrange(1)/GrdSpc_W_i(end)):floor(xrange(2)/GrdSpc_W_i(end)))*CST.DX_W/CST.DT_W,2,'omitmissing'),(1:size(CompVelWater_i.delta_x,1))*-GrdSpc_W_i(end)*CST.DX_W,'DisplayName','automated PIV','LineWidth',3)
            % plot(mean(CompVelWater.delta_x(:,floor(xrange(1)/GrdSpc_W(end)):floor(xrange(2)/GrdSpc_W(end)))*CST.DX_W/CST.DT_W,2,'omitmissing'),(1:size(CompVelWater.delta_x,1))*-GrdSpc_W(end)*CST.DX_W,'DisplayName','automated PIV','LineWidth',3)
            hold on
            %manual results
            % x = 2100 to 2600 PairNum 0432 exp5 run 2
            % u_man = [1.69e-02, 0.02624490909, 0.02624490909, 0.02624490909, 0.1181020909, 0.1106035455, 0.2005860909, 0.1612187273, 0.02811954545, 0.007498545455];
            % z_man = [-3.50E-02,-3.39E-02 -3.41E-02 -3.39E-02 -3.09E-02 -3.10E-02 -2.84E-02 -3.02E-02 -3.33E-02 -3.72E-02];

            % x = 3000 to 3500 PairNum 0259 exp5 run 2
            % u_man = [0.0	0	0.01312245455	0.01312245455	0.01312245455	0.009373181818	0.007498545455	0.009373181818	0.01124781818	0.01312245455];
            % z_man = [-4.40E-02	-3.57E-02	-3.39E-02	-3.19E-02	-3.10E-02	-3.01E-02	-2.94E-02	-2.88E-02	-2.80E-02	-3.23E-02];

            % x = 3000 to 3500 PairNum 1499 exp5 run 2
            % u_man = [9.37E-03	0.01874636364	0.005623909091	0.005623909091	0.007498545455	0.005623909091	0.01124781818	0.009373181818	0.003749272727	0.005623909091];
            % z_man = [-1.20E-01	-1.01E-01	-8.95E-02	-8.52E-02	-8.11E-02	-7.70E-02	-6.32E-02	-5.80E-02	-5.50E-02	-5.22E-02];

            % x = 1700 to 2200 PairNum 0432 exp5 run 2
            % u_man = [1.31E-02 5.62e-03 5.44E-02 1.78E-01 1.93E-01 2.12E-01 2.01E-01 9.37E-03 2.81E-02];
            % z_man = [-3.59E-02 -3.75E-02 -3.21E-02 -2.99E-02 -3.07E-02 -2.86E-02 -3.08E-02 -3.47E-02 -3.36E-02];

            % x = 3000 to 3500 PairNum 0399 exp5 run 2
            u_man = [7.50E-03 0.01312245455 0.02437027273 0.1312245455 0.1330991818 0.1349738182 0.1199767273 0.1330991818 0.1349738182 0.04874054545 0.1330991818];
            z_man = [-0.035736193 -0.03464328 -0.032766769 -0.030147902 -0.029405546 -0.028539464 -0.029591135 -0.028786916 -0.029343683 -0.030271628 -0.029343683];
            plot(u_man,z_man,'.r','DisplayName','manual','MarkerSize',25)
            legend('Location', 'southeast','Interpreter','latex')
            s = sprintf('%s %s PairNum %s, x = %d to %d pixels (%.2f to %.2f cm)',expName(1:6),runName,PairNum,xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
            title(s,'Interpreter','latex')
            set(gca,'FontSize',20,'TickLabelInterpreter','latex');
            xlabel('u (m/s)','Interpreter','latex');
            ylabel('z (m)','Interpreter','latex');
            ylim([-0.05,-0.027])

        %% Generate .tif files to feed ML PIV
        % PIV1 = double(PIV1_W);
        % PIV2 = double(PIV2_W);
        % 
        [h, w] = size(PIV1); % Image height and width
        [X1,Y1] = meshgrid([1:w], [1:h]);

        % U1 = interp2(CompVelWater.xPIV,CompVelWater.zPIV,smoothn(CompVelWater.delta_x,'robust'),X1,Y1);
        % V1 = -interp2(CompVelWater.xPIV,CompVelWater.zPIV,smoothn(CompVelWater.delta_z,'robust'),X1,Y1);
        
        PIV1 = interp2(1:size(PIV1,2),(1:size(PIV1,1))',PIVWater1_CamAngle,X1-shearGuess/2,Y1,'*linear');
        PIV2 = interp2(1:size(PIV2,2),(1:size(PIV2,1))',PIVWater2_CamAngle,X1+shearGuess/2,Y1,'*linear');
        % 
        % PIV1 = interp2(1:size(PIV1,2),(1:size(PIV1,1))',PIVWater1_CamAngle,X1-U1/2,Y1-V1/2,'*linear');
        % PIV2 = interp2(1:size(PIV2,2),(1:size(PIV2,1))',PIVWater2_CamAngle,X1+U1/2,Y1+V1/2,'*linear');
        
        % PIV1 = interp2(1:size(PIV1,2),(1:size(PIV1,1))',PIV1,X1-U1/2,Y1 - V1/2,'*linear');
        % PIV2 = interp2(1:size(PIV2,2),(1:size(PIV2,1))',PIV2,X1+U1/2,Y1 + V1/2,'*linear');
        % th = 120;
        % PIV1(PIV1 < th) = 0;
        % PIV2(PIV2*0.85 < th) = 0;

        % PIV1 = PIVWater1_CamAngle;
        % PIV2 = PIVWater2_CamAngle;

        xrange = [3000,3500];
        ysurf = 665; %677 for 0399, 665 for 0432, %665 for 0259, %800 for 1499
        xl = xrange;
        shiftmps = 0;
        shift = uint16(round(shiftmps*CST.DT_W/CST.DX_W));
        % imagesc(uint8(PIV1_W*0.84),[0,255]);
        % imagesc(uint8(PIV2_W*0.64),[0,255]);
        % set(gca,'XLim',[3000,3500],'YLim',[600,1100]
        
        % imwrite(uint8(PIV2_W(ysurf:2700,(xl(1)+shift):(xl(2)+shift))*0.64),[ImageNum_Water2, '.tif'])
        % imwrite(uint8(PIV1_W(ysurf:2700,xl(1):xl(2))*0.84),[ImageNum_Water1, '.tif'])
        
        imwrite(uint8(PIV2(ysurf:2000,(xl(1)+shift):(xl(2)+shift)).*Mask2_W(ysurf:2000,xl(1):xl(2))),[ImageNum_Water2, '.tif']) %0.64
        imwrite(uint8(PIV1(ysurf:2000,xl(1):xl(2)).*Mask1_W(ysurf:2000,xl(1):xl(2))),[ImageNum_Water1, '.tif']) %0.84
        %%
        outfname = ['/home/surflab/GitRepos/piv_liteflownet-pytorch/images/demo/DemoOutput/PIV-LiteFlowNet-en/-0_2/flow/', ImageNum_Water1, '_out.flo']; %For Pairnum 0399
        [u, v] = read_flo_file(outfname);
        figure(20)
        hold on
        % s = sprintf('CNN (piv-liteflownet-pytorch) with %.3f m/s shift',shiftmps);
        s = sprintf('CNN (piv-liteflownet-pytorch) with initial shear layer guess',shiftmps);
        plot(mean((u+double(shift)+shearGuess(ysurf:2000,(xl(1)+shift):(xl(2)+shift)))*CST.DX_W/CST.DT_W,2,'omitmissing'),(ysurf:(ysurf+size(u,1)-1))'*-1*CST.DX_W,'--b','DisplayName',s,'LineWidth',3)
        %% Plot manual and PIV velocity profiles on top of each other for either side of a trough
        % figure
        xrange = xl;
        sauto = sprintf('automated PIV (x = %d to %d pixels (%.2f to %.2f cm)',xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
        plot(mean(CompVelWater.delta_x(:,floor(xrange(1)/GrdSpc_W(end)):floor(xrange(2)/GrdSpc_W(end)))*CST.DX_W/CST.DT_W,2,'omitmissing'),(1:size(CompVelWater.delta_x,1))*-GrdSpc_W(end)*CST.DX_W,'-b','DisplayName',sauto,'LineWidth',3)
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
        
        % x = 1700 to 2200 PairNum 0432 exp5 run 2
        % u_man = [1.31E-02 5.62e-03 5.44E-02 1.78E-01 1.93E-01 2.12E-01 2.01E-01 9.37E-03 2.81E-02];
        % z_man = [-3.59E-02 -3.75E-02 -3.21E-02 -2.99E-02 -3.07E-02 -2.86E-02 -3.08E-02 -3.47E-02 -3.36E-02];


        sman = sprintf('manual (x = %d to %d pixels (%.2f to %.2f cm))',xrange(1), xrange(2), xrange(1)*CST.DX_W*1e2, xrange(2)*CST.DX_W*1e2);
        plot(u_man,z_man,'.r','DisplayName',sman,'MarkerSize',25)
        legend('Location', 'southeast','Interpreter','latex')
        s = sprintf('%s %s PairNum %s',expName(1:6),runName,PairNum);
        title(s,'Interpreter','latex')
        set(gca,'FontSize',20,'TickLabelInterpreter','latex');
        xlabel('u (m/s)','Interpreter','latex');
        ylabel('z (m)','Interpreter','latex');
        ylim([-0.05,-0.0265])
        xlim([0,0.22])
        %% warp image according to ML PIV velocity field
        tif1 = double(imread([ImageNum_Water1, '.tif']));
        tif2 = double(imread([ImageNum_Water2, '.tif']));
        [H,L] = size(tif1);
        [XT,YT] = meshgrid(1:L,1:H);

        tif1_D = interp2(XT,YT,tif1,XT-u/2,YT-v/2);
        tif2_D = interp2(XT,YT,tif2,XT+u/2, YT+v/2);
        
        
%%        


