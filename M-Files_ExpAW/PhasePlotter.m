%% Define Paths
clear
clc
% close all

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
WaterPath = [ResultsPath 'Water/'];
SaveSurfWaterPath = [WaterPath 'Surfaces/'];

runName = sprintf('Run%d',runNum);

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water
%% Prep for surface detection and everything that relies on that.
frames = 0:1:1499; %numbers of frames to loop through. First frame in the experiment is denoted by index 0. Must end on odd number so that you get an integer number of pairs
PairNumsInt = (ceil(frames(1)/2)):(floor(frames(end)/2));
numPairs = length(PairNumsInt);

nF = length(frames);
fXs = zeros(nF, 3639);
fYs = zeros(nF, 3639);
XPIVW_PIVSurfW_Surfaces = NaN(nF,5726);
PIVW_PIVSurfW_Surfaces = NaN(nF,5726);
surfVelMeans = NaN(numPairs,1);
shearLayerDepths = NaN(numPairs,1);

mpp = 6.493178e-5; % meters per pixel, for main dataset surface images only

pps = 14.5; %pairs per second

spp = 1/pps; % seconds per pair

yLimits = [1950,2150];
dy = yLimits(2)-yLimits(1);

dt_pair = 22.222e-3; % time between pictures in a given pair

dydt = dy/dt_pair;

DeltaT = nF/2*spp;

% FI = 0;%First index
% LI = length(PIVWaterDir)-1;
% image_index = FI+1:LI; %1, 3, 5,... Set of indices to loop through. Images are processed in pairs, hence the increment of 2

%% Detect surface in every frame specified by frames
redetectSurf = false; % set to true if you want to run surface detection here. Set to false if you want to read the surf from the Water/Surfaces .mat files

if redetectSurf
    tic
    parfor framesCtr = 1:length(frames)
        idx = frames(framesCtr);
        PIV1Dir_temp = PIVWaterDir;
    
        imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
        imagename
        
        [IM1] = load_Image_IOCoreView_12MP(imagename);
        PIVSurf_W1_Raw = IM1;
        PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';
    
        [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
        [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);
    
        PIV1_W = [];
        [BadFramePIVSurfW1, XPIVSurfW1_Surface, PIVSurfW1_Surface, XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, PIVW1_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW1_CamAngle,PIV1_W,'1');
        
    
    
        % imagesc(PIVSurfW1_CamAngle, [0,70])
        % hold on
        % set(gca,'DataAspectRatio',[1 1 1])
        % ylim([1900,2200])
        % plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
        % pause(0.1)
        
        fXs(framesCtr,:) = XPIVSurfW1_Surface;
        fYs(framesCtr,:) = PIVSurfW1_Surface;
        XPIVW_PIVSurfW_Surfaces(framesCtr,:) = XPIVW_PIVSurfW1_Surface;
        PIVW_PIVSurfW_Surfaces(framesCtr,:) = PIVW_PIVSurfW1_Surface;
    end
    toc
    save('temp.mat','-v7.3')
    clear
    load('temp.mat')
else
    tic
    for i = 1:numPairs
        PairNumInt = PairNumsInt(i);
        PairNum = sprintf('%04d',PairNumInt)
        SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];

        try
            load([SaveSurfWaterPath SurfFileName '.mat']);

            fXs(i*2-1,:) = SurfRes_Water1.XPIVSurfW1_Surface;
            fYs(i*2-1,:) = SurfRes_Water1.PIVSurfW1_Surface;
            XPIVW_PIVSurfW_Surfaces(i*2-1,:) = PixRes_Water1.XPIVW_PIVSurfW1_Surface;
            PIVW_PIVSurfW_Surfaces(i*2-1,:) = PixRes_Water1.PIVW_PIVSurfW1_Surface;
            fXs(i*2,:) = SurfRes_Water2.XPIVSurfW2_Surface;
            fYs(i*2,:) = SurfRes_Water2.PIVSurfW2_Surface;
            XPIVW_PIVSurfW_Surfaces(i*2,:) = PixRes_Water2.XPIVW_PIVSurfW2_Surface;
            PIVW_PIVSurfW_Surfaces(i*2,:) = PixRes_Water2.PIVW_PIVSurfW2_Surface;
            shearLayerDepths(i) = (SurfVel.bottomShearLayer-SurfVel.mwl)*CST.DX_W;
            surfVelMeans(i) = mean(SurfVel.delta_x,'all','omitmissing')*CST.DX_W/CST.DT_W;

        catch
            SurfFileName = [expName '_Scene' runName '_Surfaces_' PairNum];
            disp(['Could not load surface from file ' SaveSurfWaterPath SurfFileName '.mat'])
        end
    end
    toc
end
%% Assemble composite image for hovmoller plot
tic 
CompImg = NaN(ceil(DeltaT*dydt),4176);
for framesCtr = 1:length(frames)
    idx = frames(framesCtr);

    PIV1Dir_temp = PIVWaterDir;

    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    imagename
    
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';

    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    [PIVSurfW1_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);

    if mod(idx,2) == 0
        CompImg(round((framesCtr-1)/2*spp*dydt)+1:round((framesCtr-1)/2*spp*dydt)+dy+1,1:size(PIVSurfW1_CamAngle,2)) = PIVSurfW1_CamAngle(yLimits(1):yLimits(2),:)/mean(PIVSurfW1_CamAngle,'all');
    else
        CompImg(round(((framesCtr-2)/2*spp+dt_pair)*dydt)+1:round(((framesCtr-2)/2*spp+dt_pair)*dydt)+dy+1,1:size(PIVSurfW1_CamAngle,2)) = PIVSurfW1_CamAngle(yLimits(1):yLimits(2),:)/mean(PIVSurfW1_CamAngle,'all');
    end
end
toc

%% FULL SYNOPSIS: Setup for plotting eta_var, eta_x_var, windspeed, surface speed, Part 1: Calculate eta stats
trange = [20, 37];
PairNumRange = trange/spp;
t = zeros(1,nF);
t(1) = floor(frames(1)/2)*spp + mod(frames(1),2)*dt_pair;

for i = 2:nF
    if mod(frames(i),2)==1
        t(i) = t(i-1)+dt_pair;
    else
        t(i) = t(i-1)+spp-dt_pair;
    end
end
PairNumCont = t/spp;
CST = load('CST.mat'); %This is bad and should be changed because it relies on a file in the main M-Files_ExpAW directory that is specific to a give experiment
etaMeanPIVW = mean(PIVW_PIVSurfW_Surfaces*CST.DX_W, 'all','omitmissing');

eta = (fYs-mean(fYs(1:20,:),'all'))*mpp;
eta = detrend(eta);
% eta = highpass(eta',4,1/mpp,Steepness=0.999999)';
x_eta = fXs*mpp;
eta_var = sum((eta-mean(eta,2)).^2,2)/(size(eta,2)-1);

eta_x = diff(eta/mpp,1,2);

eta_x_var = sum((eta_x-mean(eta_x,2)).^2,2)/(size(eta_x,2)-1);

%% FIND PHASE SPEED: using cross-correlation method

phaseSpeed = nan(numPairs,1);
P_threshold = 0.0004;

for i = 1:numPairs
    fa = (i-1)*2 + 1;
    fb = (i-1)*2 + 2;
    etaa = eta(fa,:);
    etab = eta(fb,:);
    [~,P] = islocalmin(etaa);
    if max(P) > P_threshold
        [r,lags] = xcorr(etaa,etab);
        [~,iMaxR] = max(r);
        lag = lags(iMaxR)
        PS = -lag*mpp/dt_pair
        phaseSpeed(i) = PS;
        % figure(4)
        % hold off
        % plot(((lag+1):(length(etab) + lag))*mpp, -etab)
        % hold on
        % daspect([1,1,1])
        % plot((1:length(etab))*mpp,-etaa)
        % ylim([-5e-3,5e-3])
        % pause
    end

end

%% MAKE VIDEO FRAMES with phase speed label (rewritten with Claude)
clear frames
% Determine which indices we're processing (even only)
allIdx = 350:(size(fXs,1)-400);
evenIdx = allIdx(mod(allIdx, 2) == 0);
numFrames = length(evenIdx);

% Pre-allocate cell array for frames
frames = cell(numFrames, 1);

% Figure dimensions matching your current setup
figWidth = 2000;
figHeight = 300;

parfor k = 1:numFrames
    idx = evenIdx(k);
    pairNum = floor(idx/2);
    
    % Create figure with specific pixel dimensions
    f1 = figure('Units', 'pixels', 'Position', [0, 0, figWidth, figHeight], ...
                'Color', 'w', 'Visible', 'off');  % 'off' can speed things up
    
    % Make axes fill the entire figure (no border)
    ax = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    
    imagename = [PIVWaterDir(idx+1).folder '/' PIVWaterDir(idx+1).name];
    
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    PIVSurf_W1_Raw = IM1;
    PIVSURF_W1 = PIVSurf_W1_Raw./(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)/max(smooth(mean(PIVSurf_W1_Raw(2000:end,:)),1000)))';
    
    [PIVSurfW1_Undistorted] = PIVSurfW_LensDistCorr(PIVSURF_W1);
    [PIVSurf_CamAngle] = PIVSurfWater_CamAngle_Correction(PIVSurfW1_Undistorted);
    surfImg = PIVSurf_CamAngle;
    
    imagesc(ax, surfImg, [0, 60])
    hold(ax, 'on')
    axis(ax, 'off')
    daspect(ax, [1, 1, 1])
    colormap(ax, gray)
    plot(ax, fXs(idx+1, 1:3:end), fYs(idx+1, 1:3:end), '-r', 'LineWidth', 2)
    
    xl = [501, 4139];
    yl = [1750, 2250];
    
    xlim(ax, xl);
    ylim(ax, yl);
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm / mpp;
    xsb = [xl(1) + (xl(2)-xl(1))*0.05, xl(1) + (xl(2)-xl(1))*0.05 + lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1) * [1 1];
    
    plot(ax, xsb, ysb, '-k', 'LineWidth', 10);
    
    sbl = sprintf('%d cm', lsbm*100);
    text(ax, xsb(2) + (xl(2)-xl(1))*0.01, ysb(2), sbl, 'FontSize', 34, 'Interpreter', 'latex');
    
    text(ax, 0.85*xsb(1), (yl(2) - (yl(2)-yl(1))*0.3), 'Water', 'FontSize', 34, 'Interpreter', 'latex')
    text(ax, 0.85*xsb(1), (yl(2) - (yl(2)-yl(1))*0.6), 'Air', 'FontSize', 34, 'Interpreter', 'latex', 'Color', [1,1,1])
    text(ax, 0.85*xsb(1), (yl(2) - (yl(2)-yl(1))*0.8), '$Wind \rightarrow$', 'FontSize', 34, 'Interpreter', 'latex', 'Color', [1,1,1])
    
    ttext = sprintf('t = %.1f s', t(idx+1)-11);
    text(ax, xl(1) + (xl(2)-xl(1))*0.92, yl(2) - (yl(2)-yl(1))*0.9, ttext, 'FontSize', 34, 'Interpreter', 'latex', 'Color', [1,1,1])
    
    if ~isnan(phaseSpeed(pairNum+1))
        cptext = sprintf('$c_p = %.1f$ cm/s', phaseSpeed(pairNum+1)*100);
        text(ax, xl(1) + (xl(2)-xl(1))*0.45, yl(2) - (yl(2)-yl(1))*0.9, cptext, 'FontSize', 34, 'Interpreter', 'latex', 'Color', [1,1,1])
    end
    
    drawnow
    
    % Capture frame to cell array instead of saving to disk
    frame = getframe(f1);
    frames{k} = frame.cdata;  % Store just the RGB image data
    
    close(f1)
end

%% Write video sequentially (this is fast compared to frame generation) with Claude help
fprintf('Writing video...\n');

vidObj = VideoWriter(['videoframes/' DataPath(end-23:end-1) '.avi'], 'Uncompressed AVI');
vidObj.FrameRate = 14.5/4;
open(vidObj);

for k = 1:numFrames
    writeVideo(vidObj, frames{k});
end

close(vidObj);
aviFile = ['videoframes/' DataPath(end-23:end-1) '.avi'];
mp4File = ['videoframes/' DataPath(end-23:end-1) '.mp4'];
system(['ffmpeg -i "' aviFile '" -c:v libx264 -pix_fmt yuv420p -crf 18 "' mp4File '"']);
fprintf('Done!\n');



% Clear frames from memory
% clear frames;

%% FULL SYNOPSIS: Setup for plotting eta_var, eta_x_var, windspeed, surface speed, Part 2: Calculate Air PIV stats
withAir = true;

PairNums = floor(t(1)/spp):floor(t(end)/spp);

if withAir
    UAirDir = [ResultsPath 'Air/CALCULATED_FIELDS/Cartesian Fields/Velocity/'];
    
    PairNum = PairNums(281);
    fname = sprintf('%s_CartesianAir_%04d.mat',expRunName(1:end-1), PairNum);
    Cartesian_Air = load([UAirDir,fname]);
    CST = Cartesian_Air.CST;
    
    up_b = NaN(length(PairNums),size(Cartesian_Air.Mask,1));
    wp_b = NaN(length(PairNums),size(Cartesian_Air.Mask,1));
    u_mean = NaN(length(PairNums),1);
    w_mean = NaN(length(PairNums),1);
    upwp_b = NaN(length(PairNums),size(Cartesian_Air.Mask,1));
    dupwpdz_b = NaN(length(PairNums),size(Cartesian_Air.Mask,1)-1);
    
    tic
    parfor n = 1:length(PairNums)
        try
            PairNum = PairNums(n);
            fname = sprintf('%s_CartesianAir_%04d.mat',expRunName(1:end-1), PairNum);
        
            Cartesian_Air = load([UAirDir,fname]);
            CST = Cartesian_Air.CST;
            PairNum = str2double(Cartesian_Air.PairNum);
        
            u = Cartesian_Air.u.*Cartesian_Air.Mask*CST.DX/CST.DT;
            w = Cartesian_Air.w.*Cartesian_Air.Mask*CST.DX/CST.DT; %Positive is up for w and z. Postive is down for v and y.
            speed = (u.^2 + w.^2).^0.5;
            mask = Cartesian_Air.Mask;
            
            u_b = mean(u,2,'omitnan');
            w_b = mean(w,2,'omitnan');
            
            up = u - u_b;
            wp = w - w_b;
        
            up_b(n,:) = mean(up,2,'omitnan');
            wp_b(n,:) = mean(wp,2,'omitnan');
            
            upwp_b(n,:) = mean(up.*wp,2,'omitnan');
            
            dupwpdz_b(n,:) = mean(diff(up.*wp,1,1),2,'omitnan');
    
            u_mean(n) = mean(u,'all','omitmissing');
            w_mean(n) = mean(w,'all','omitmissing');
        
            disp(['Pair ' Cartesian_Air.PairNum ' Finished!'])
        catch
            disp('Missing Pair')
    
        end
    
    end
    toc
    
    etaMeanPIVA = mean(Cartesian_Air.Surface,'all','omitmissing')*CST.DX;
    
    
    plot(PairNums, movmean(upwp_b(:,floor((etaMeanPIVA-0.01)/CST.DX)),10))
end

t_airPIV = PairNums*spp;

%% FULL SYNOPSIS: Plot eta_var, eta_x_var, windspeed, surface speed
withMan = false;

co = colororder('gem');
co(6,:) = [0.15 0.15 0.8];
co(5,:) = [0.1 0.4 0.3];

% Tile 1, eta_var
% Start and end times for the initial growth phase exponential fit
if runNum == 2 || runNum == 3
    switch ExpAW
        case '1'
            tfit_i = 43;
            tfit_f = 45.7;
        case '2'
            tfit_i = 33;
            tfit_f = 36.6;
        case '3'
            tfit_i = 28;
            tfit_f = 33;
        case '4'
            tfit_i = 32.5;
            tfit_f = 35.5;
        case '5'
            tfit_i = 27;
            tfit_f = 29.58;
        case '7'
            tfit_i = 26.2;
            tfit_f = 30.3;
        case '8'
            tfit_i = 28;
            tfit_f = 32;
    end
end

% tfit_i = 32.5;
% tfit_f = 35.2;

ir = [0,0];
[~,i] = min(abs(t-tfit_i));
ir(1) = i;
[~,i] = min(abs(t-tfit_f)); 
ir(2) = i;

eta_x_var(eta_var > 5e-6) = nan;
eta_x_var = smoothn(eta_x_var, 1e-9,'robust');

eta_var(eta_var > 5e-6) = nan;
eta_var = smoothn(eta_var, 1e-9,'robust');
eta_var_fit = fit(t(ir(1):ir(2))',eta_var(ir(1):ir(2)),'exp1');

figure(4)
if withAir
    numTiles = 4;
else
    numTiles = 3;
end
    
tlay = tiledlayout(numTiles,1);

tileNum = 1;

ax1_1 = axes(tlay);
ax1_1.Layout.Tile = tileNum;
ax1_1.XAxisLocation = 'bottom';

p1_1 = plot(ax1_1,t,eta_var,'LineWidth',3,'Color',co(1,:));
hold on
ax1_1.YColor = co(1,:);
% xlabel(ax1_1,'Time (s)','Interpreter','latex')
ylabel('$\mathrm{Var}[\eta]\ \mathrm{(m^2)}$','Interpreter','latex')
set(ax1_1,'FontSize',18)
set(ax1_1,'TickLabelInterpreter','latex')
set(ax1_1, 'YLim',[0,2.5e-7],'XLim',trange)
ax1_1.Box = 'off';
s = [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)];
title(s,'Interpreter','latex')

s = sprintf('Points %.2f s to %.2f s',tfit_i,tfit_f);
p1_3 = plot(ax1_1,t(ir(1):ir(2))',eta_var(ir(1):ir(2)),'.k', 'MarkerSize',15,'DisplayName',s);

fit_line_t = (t(ir(1)):0.01:t(ir(2)))';
fit_line = eta_var_fit.a.*exp(fit_line_t.*eta_var_fit.b);
s = sprintf('Fit %.2f s to %.2f s, Growth Rate %.2f 1/s',tfit_i,tfit_f,eta_var_fit.b);
p1_4 = plot(ax1_1,fit_line_t, fit_line, '-r', 'LineWidth',3,'DisplayName',s);
hold off

% Tile 1, eta_x_var
ax1_2 = axes(tlay);
hold on
ax1_2.YColor = co(2,:);
ax1_2.Color = 'none';
ax1_2.Box = 'off';
ax1_2.Layout.Tile = tileNum;
ax1_2.YAxisLocation = 'right';
p1_2 = plot(ax1_2,t,eta_x_var,'LineWidth',3,'Color',co(2,:));
set(ax1_2,'FontSize',18)
tickint = 10;
xticks(ax1_2,t(1:tickint:end))
xticklabels(ax1_2,PairNumCont(1:tickint:end))
ylabel(ax1_2,'$\mathrm{Var}[\eta_x]$','Interpreter','latex')
ax1_2.XAxisLocation = 'top';
ax1_2.XAxis.FontSize = 6;
ax1_2.XAxis.Color = [0.5,0.5,0.5];
set(ax1_2,'TickLabelInterpreter','latex')
set(ax1_2, 'YLim',[0,0.05],'XLim',trange)
tempyl = ylim;
xlabel(ax1_2,'PairNum','Position',[trange(2) - (trange(2)-trange(1))*0.015,tempyl(2)*1.07],'Interpreter','latex')
% legend([p1,p2])
legend(ax1_1,[p1_3,p1_4],'Interpreter','latex','Location','northwest')

hold off

if withAir
    u_filtWndw = 50;
    %Tile 2: Wind Speed
    tileNum = tileNum + 1;
    ax2_1 = axes(tlay);
    hold on
    ax2_1.Layout.Tile = tileNum;
    ax2_1.Box = 'off';
    ax2_1.YColor = co(6,:);
    ylabel('$\mathrm{Var}[\eta]\ \mathrm{(m^2)}$','Interpreter','latex')
    set(ax1_1,'FontSize',18)
    set(ax1_1,'TickLabelInterpreter','latex')
    set(ax1_1, 'YLim',[0,2.5e-7],'XLim',trange)
    p2_1 = plot(ax2_1,t_airPIV, movmean(u_mean,1,'omitmissing'),'-','LineWidth',3,'Color',co(6,:),'DisplayName','$\overline{U}_{air}\ \mathrm{(m/s)}$');
    ylabel(ax2_1,'$\overline{U}_{air}\ \mathrm{(m/s)}$','Interpreter','latex')
    set(ax2_1,'FontSize',18)
    set(ax2_1,'TickLabelInterpreter','latex')
    set(ax2_1, 'YLim',[0,6],'XLim',trange)
    ax2_1.Box = 'off';
    


    u_mean_filt = movmean(u_mean,u_filtWndw,'omitmissing');
    [umax, i_umax] = max(u_mean_filt);
    umax
    % [~,i_acc_f] = min(abs(u_mean_filt-0.9*umax));
    % [~,i_acc_i] = min(abs(u_mean_filt-0.05*umax));
    
    i_acc_f = find(u_mean_filt > 0.9*umax);
    i_acc_f = i_acc_f(1);
    i_acc_i = find(u_mean_filt > 0.05*umax);
    i_acc_i = i_acc_i(1);

    t_airfit = t_airPIV(i_acc_i:i_acc_f);
    u_airfit = u_mean(i_acc_i:i_acc_f);
    t_airfit = t_airfit(isfinite(u_airfit));
    u_airfit = u_airfit(isfinite(u_airfit));
    airVelFit = polyfit(t_airfit,u_airfit',1)
    air_acc = airVelFit(1);
    s = sprintf("Linear Fit ($d\\overline{U}_{air}/dt = %.2f\\ \\mathrm{m/s^2}$), $t = %.2f$ to $%.2f$s",air_acc, t_airfit(1),t_airfit(end));
    p2_4 = plot(ax2_1,t_airfit,polyval(airVelFit,t_airfit),'--','LineWidth',3,'Color',co(6,:),'DisplayName',s');
    legend(p2_4,'Interpreter', 'latex','Location', 'northwest')
    hold off

    % figure(5)   
    % plot(t_airPIV(2:end),movmean(diff(movmean(u_mean,50,'omitmissing')),50)./diff(t_airPIV))
    % hold on
    % ylabel('Acceleration of $\overline{U}_{air}\ \mathrm{(m/s^2)}$','Interpreter','latex')
    % xlabel('Time (s)')
    % hold off

    figure(4)
    
    ax2_2 = axes(tlay);
    ax2_2.Layout.Tile = tileNum;
    ax2_2.Color = 'none';
    ax2_2.YAxis.Visible = 'off';
    ax2_2.Box = 'off';
    hold on
    xticks(ax2_2,t(1:tickint:end))
    xticklabels(ax2_2,PairNumCont(1:tickint:end))
    set(ax2_2,'XLim',trange)
    ax2_2.XAxisLocation = 'top';
    ax2_2.XAxis.FontSize = 6;
    ax2_2.XAxis.Color = [0.5,0.5,0.5];
    set(ax2_2,'TickLabelInterpreter','latex')
    set(ax2_2, 'YLim',[0,0.05],'XLim',trange)
    hold off
    
    upwp_b_plot_loc = 0.01; %height above mwl (m)
    ax2_3 = axes(tlay);
    ax2_3.Layout.Tile = tileNum;
    hold on
    ax2_3.Color = 'none';
    ax2_3.YAxis.Visible = 'on';
    ax2_3.XAxis.Visible = 'off';
    ax2_3.Box = 'off';
    ax2_3.YAxisLocation = 'right';
    ax2_3.YColor = co(7,:);
    p2_3 = plot(ax2_3,t_airPIV, 100*sqrt(-movmean(upwp_b(:,floor((etaMeanPIVA-upwp_b_plot_loc)/CST.DX)),10,'omitmissing')),'LineWidth',3,'Color',co(7,:),'DisplayName',"$\sqrt{-\overline{u'w'}}$");
    ylabel(ax2_3,"$\sqrt{-\overline{u'w'}}\ \mathrm{(cm/s)}$",'Interpreter','latex')
    set(ax2_3,'FontSize',18)
    set(ax2_3,'TickLabelInterpreter','latex')
    set(ax2_3, 'XLim',trange)
    hold off
end

% Tile 3: Surface Velocity
tileNum = tileNum + 1;
if withMan
    manStats = readtable([ManualResultsPath, expRunName(1:end-1), '_ManualDataSummary.csv']);
    manPairNums = manStats{:,1};
    manTimes = manPairNums*spp;
    manSurfVels = manStats{:,2};
    manShearDepths = manStats{:,3}+etaMeanPIVW;
end

[maxSurfVel, ImaxSurfVel] = max(movmean(surfVelMeans,5));
maxSurfVel
disp("t at maxsurfvel (s): " + t_airPIV(ImaxSurfVel))
[linStartSurfVel, IlinStartSurfVel] = min(abs(surfVelMeans(1:ImaxSurfVel) - maxSurfVel*0.1));

surfVelMeansFitMask = false(size(surfVelMeans));
surfVelMeansFitMask(IlinStartSurfVel:ImaxSurfVel) = true;
surfVelMeansFitMask = surfVelMeansFitMask & isfinite(surfVelMeans);
t_lf = t_airPIV(surfVelMeansFitMask);
surfVelFit = polyfit(t_lf,surfVelMeans(surfVelMeansFitMask),1);
A = surfVelFit(1);
ySurfVelFit = polyval(surfVelFit,t_lf);
D_msv98 = NaN(length(t_lf),1); % shear layer depth given by Melville,Shear,Veron 1998 for a linearly increasing surface velocity
t_0 = -surfVelFit(2)/A;
for n = 1:length(t_lf)
    t_rel = t_lf(n)-t_0;
    syms z
    et = -z/(2*(CST.WATER_DVISCOSITY/CST.WATER_DENSITY*t_rel)^0.5);
    D_msv98(n) = vpasolve(0.1 == ((1+2*et^2)*erfc(et)-2/(pi^0.5)*et*exp(-et^2)),z);
end

ax3_1 = axes(tlay);
hold on
ax3_1.Layout.Tile = tileNum;
ax3_1.Box = 'off';
p3_2 = plot(ax3_1,t_airPIV,surfVelMeans*100,'-','LineWidth',3,'Color',co(4,:),'DisplayName','Automated $U_{surf}$');
s = sprintf("Linear Fit ($dU_{surf}/dt = %.2f\\ \\mathrm{cm/s^2}$)",100*surfVelFit(1));
p3_3 = plot(ax3_1,t_lf,ySurfVelFit*100,'--','LineWidth',3,'Color',co(4,:),'DisplayName',s);
if withMan
    p3_1 = plot(ax3_1,manTimes,manSurfVels*100,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',co(4,:),'DisplayName','Manual $U_{surf}$');
end
p3_7 = plot(ax3_1,t_airPIV, phaseSpeed*100,':','LineWidth',3,'Color',co(4,:),'DisplayName','Cross-correlation $c_p$')
ylim([0,35])
set(ax3_1,'FontSize',18)
set(ax3_1,'TickLabelInterpreter','latex')
set(ax3_1,'XLim',trange)
ylabel(ax3_1,'$c_p,\ U_{surf}$ (cm/s)','Interpreter','latex')
ax3_1.YColor = co(4,:);

ax3_2 = axes(tlay);
ax3_2.Layout.Tile = tileNum;
ax3_2.Color = 'none';
ax3_2.YAxis.Visible = 'off';
hold on
xticks(ax3_2,t(1:tickint:end))
xticklabels(ax3_2,PairNumCont(1:tickint:end))
set(ax3_2,'XLim',trange)
ax3_2.XAxisLocation = 'top';
ax3_2.XAxis.FontSize = 6;
ax3_2.XAxis.Color = [0.5,0.5,0.5];
set(ax3_2,'TickLabelInterpreter','latex')
hold off

ax3_3 = axes(tlay);
ax3_3.Layout.Tile = tileNum;
hold on
ax3_3.Color = 'none';
ax3_3.Box = 'off';
ax3_3.YAxisLocation = 'right';
ax3_3.XAxis.Visible = 'off';
p3_5 = plot(ax3_3,t_airPIV,-shearLayerDepths*100,'-','LineWidth',3,'Color',co(5,:),'DisplayName','Automated $D_{SL}$');
p3_6 = plot(ax3_3,t_lf,D_msv98*100,'--','LineWidth',3,'Color',co(5,:),'DisplayName','Melville, Shear, Veron (1998) $D_{SL}$');
if withMan
    p3_4 = plot(ax3_3,manTimes,manShearDepths*100,"^",'MarkerSize',10,'Color',co(5,:),'MarkerEdgeColor','k','MarkerFaceColor',co(5,:),'DisplayName','Manual $D_{SL}$');
end
disp("Shear Layer Depth at max surf vel: " + shearLayerDepths(ImaxSurfVel)*100 + " cm")
ylabel(ax3_3,{'Shear Layer Depth'; '(cm)'},'Interpreter','latex')
ax3_3.YColor = co(5,:);
set(ax3_3,'FontSize',18,'TickLabelInterpreter','latex','XLim',trange,'YLim',[-2.3,0])
if withMan
    legend([p3_1, p3_2, p3_3, p3_4,p3_5,p3_6,p3_7],'Interpreter','latex','Location', 'northwest')
else
    legend([p3_2, p3_3,p3_5,p3_6,p3_7],'Interpreter','latex','Location', 'northwest')
end

% Tile 4: Surface Area
tileNum = tileNum + 1;
ax4_1 = axes(tlay);
hold on
ax4_1.Layout.Tile = tileNum;
ax4_1.Box = 'off';
A_surf = zeros(1,length(t));
for n = 1:length(t)
    for i = 1:(size(eta,2)-1)
        A_surf(n) = A_surf(n) + ((eta(n,i+1) - eta(n,i))^2 + (x_eta(n,i+1) - x_eta(n,i))^2).^0.5;
    end
end

l_interface = x_eta(1,end) - x_eta(1,1);
p4_1 = plot(ax4_1,t,(A_surf - l_interface)/l_interface*100,'LineWidth',3,'Color','k','DisplayName', 'surface area')
hold on
xlabel(ax4_1,'Time (s)','Interpreter','latex')
ylabel(ax4_1,{'$\mathrm{Surface\ Area}$'; '$\mathrm{Increase\ (\%)}$'},'Interpreter','latex')
set(ax4_1,'FontSize',18)
set(ax4_1,'TickLabelInterpreter','latex')
set(ax4_1, 'YLim',[0,3],'XLim',trange)

ax4_2 = axes(tlay);
ax4_2.Layout.Tile = tileNum;
ax4_2.Color = 'none';
ax4_2.YAxis.Visible = 'off';
hold on
xticks(ax4_2,t(1:tickint:end))
xticklabels(ax4_2,PairNumCont(1:tickint:end))
set(ax4_2,'XLim',trange)
ax4_2.XAxisLocation = 'top';
ax4_2.XAxis.FontSize = 6;
ax4_2.XAxis.Color = [0.5,0.5,0.5];
set(ax4_2,'TickLabelInterpreter','latex')
set(ax4_2, 'YLim',[0,0.05],'XLim',trange)

if withAir
    linkaxes([ax1_1,ax1_2,ax2_1,ax2_2,ax2_3, ax3_1, ax3_2, ax3_3, ax4_1, ax4_2],'x')
else
    linkaxes([ax1_1,ax1_2, ax3_1, ax3_2, ax3_3, ax4_1, ax4_2],'x')
end

%% Plot Wavenumber Spectrum
nX = size(eta,2);
ppm = 1/mpp;

kspec_mag = abs(fft(eta-mean(eta,2),nX,2));

%Calculate contribution to total variance by each wavenumber, normalized by wavenumber bandwidth
%Sum of kspec_mag_n*(k step size) = kspec_mag_n*(2*pi/nX/mpp) is equal to
%the total variance
kspec_vareta_n = kspec_mag.^2/(nX-1)*mpp/(2*pi); 

figure(13)
hold off
imagesc(kspec_mag.^2/(nX-1)*mpp/(2*pi),'XData',2*pi*(ppm/nX*(0:nX-1)),'YData',t,[0,1e-9])
hold on
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0,500])
ylim([22,30])
xlabel('k ($\mathrm{m^{-1}}$)','Interpreter','latex')
ylabel('time (s)','Interpreter','latex')
colormap gray
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = "Wavenumber contribution to Var[$\eta$] ($\mathrm{m^3}$)";
c.Label.FontSize = 24;
c.TickLabelInterpreter = 'latex';
s = [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)];
title(s,'Interpreter','latex')
%% Plot dispersion relation for a certain time interval
T_man = (22.8843-22.3671+23.1602-22.6084)/2;
g = 9.81;
sig_man = 2*pi/T_man;
dt_man = 23.1602-22.8843;
dx_man = 0.267194-0.0339593;
lambda_man = dx_man*T_man/dt_man;


is = 328*2-50:328*2+50;%385*2:450*2;
sigma = 1:0.1:30;
k = 1:50;
x = x_eta(1,:);
eta_dm = eta - mean(eta,2,'omitmissing');
sig_spec = nufft(eta_dm(is,:),t(is),sigma/(2*pi),1);
dr = nufft(nufft(eta_dm(is,:),t(is),sigma/(2*pi),1),x,k/(2*pi),2);
figure(17)
imagesc(k, sigma, abs(dr))
hold on
plot(k,sqrt(g*k),'-r','LineWidth',3,'DisplayName','Gravity wave dispersion relation');
plot(2*pi/lambda_man,sig_man,'^','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',25,'DisplayName','Manual estimate of frequency and wavenumber')
legend('Interpreter','latex')
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
xlabel('k (rad/m)','Interpreter','latex');
ylabel('$\omega$ (rad/s)', 'Interpreter','latex');
c = colorbar;
colormap gray;
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.Label.String = '2D FFT Intensity (m)';
st = sprintf('$\\eta$ from $t =$ %.4f s to %.4f s, spatially de-meaned at each time step, for %s',t(is(1)),t(is(end)),[expRunName(1:6),' ',expRunName(end-4:end-1)]);
title(st,'Interpreter','latex');
hold off
figure(18)
imagesc(x,sigma,real(sig_spec))
figure(19)
imagesc(x,t(is),eta(is,:))
hold on
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
st = sprintf('$\\eta$ from $t =$ %.4f s to %.4f s for %s',t(is(1)),t(is(end)),[expRunName(1:6),' ',expRunName(end-4:end-1)]);
title(st,'Interpreter','latex')
xlabel('x (m)','Interpreter','latex');
ylabel('t (s)','Interpreter','latex');
c = colorbar;
colormap gray;
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
c.Label.String = '$\eta$ (m)';
hold off;
%% Plot initial growth rate spectrum

growthrates = zeros(1,nX);
parfor i = 1:nX
    sfit = fit(t(ir(1):ir(2))',kspec_vareta_n(ir(1):ir(2),i),'exp1');
    growthrates(i) = sfit.b;
end

figure(14)
% hold off
s = sprintf("%.2f s to %.2f s",tfit_i,tfit_f);
plot(2*pi*(ppm/nX*(0:nX-1)),growthrates,'LineWidth',4,'DisplayName',s)
hold on
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0,500])
xlabel('k ($\mathrm{m^{-1}}$)','Interpreter','latex')
ylabel('growth rate of wavenumber contribution to Var[$\eta$] (1/s)','Interpreter','latex')

s = [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)];
title(s,'Interpreter','latex')
legend('Interpreter','latex')

%% Determine Asymmetry of troughs
P_threshold = 7;
lambda_max = 4e-2; % asymmetry detection wavelength maximum.
slopePThreshold = 0.1;

slopeDiffMean = nan(1,length(t));
slopeDiffMedian = nan(1, length(t));
slopeDiffMode = nan(1,length(t)); 
numTroughs = nan(1,length(t));
for i = 1:length(t)
    [TF, P] = islocalmin(-fYs(i,:));

    % slopeMask = islocalmax(abs(diff(fYs(i,:))));
    % steepPtsMask(1,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(1);
    % steepPtsMask(2,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(2);
    % steepPtsMask(3,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(3);
    
    minsMask = P > P_threshold;
    imins = find(minsMask);
    xmins = fXs(i,minsMask);
    ymins = fYs(i,minsMask);


    if ~isempty(ymins)
        slope = -diff(fYs(i,:));
        [slopeMins,slopeMinsP] = islocalmin(slope);
        [slopeMaxs,slopeMaxsP] = islocalmax(slope);
    
        if mod(i,2) == 1
            offset = round((i-1)/2*spp*dydt)+1-yLimits(1);
            offset2 = round((i-1)/2*spp*dydt)+1;
        else
            offset = round(((i-2)/2*spp+dt_pair)*dydt)+1-yLimits(1);
            offset2 = round(((i-2)/2*spp+dt_pair)*dydt)+1;
        end

        slopeDiffs = nan(1,length(ymins));
        for j = 1:length(ymins)
            foundSlopeMin = false;
            foundSlopeMax = false;
            k = imins(j);
            troughSlopeMin = 0;
            troughSlopeMax = 0;
            while ~foundSlopeMin && k > 0 && imins(j)-k < lambda_max/2/mpp
                if slopeMinsP(k) > slopePThreshold
                    foundSlopeMin = true;
                    troughSlopeMin = slope(k);
                    plot(fXs(i,k)*mpp,(fYs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                end
                k = k-1;
            end
    
            k = imins(j);
            while ~foundSlopeMax && k < length(slope) && k-imins(j) < lambda_max/2/mpp
                if slopeMaxsP(k) > slopePThreshold
                    foundSlopeMax = true;
                    troughSlopeMax = slope(k);
                    plot(fXs(i,k)*mpp,(fYs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                end
                k = k+1;
            end
    
            slopeDiff = abs(troughSlopeMax + troughSlopeMin);
            slopeDiffText = sprintf('SlopeDiff %.2f',slopeDiff);
            slopeDiffs(j) = slopeDiff;
        end
        slopeDiffMean(i) = mean(slopeDiffs);
        slopeDiffMedian(i) = median(slopeDiffs);
        slopeDiffMode(i) = mode(slopeDiffs);
        numTroughs(i) = length(ymins);
    end
end


figure(4)
ax4_3 = axes(tlay);
ax4_3.Layout.Tile = tileNum;
hold on
ax4_3.Color = 'none';
ax4_3.Box = 'off';
ax4_3.YAxisLocation = 'right';
ax4_3.XAxis.Visible = 'off';
p4_2 = plot(t,slopeDiffMean,'LineWidth',2,'Color',co(1,:),'DisplayName', '$$\Delta \eta_x$$ mean');
p4_3 = plot(t,slopeDiffMedian,'LineWidth',2,'Color',co(2,:),'DisplayName', '$\Delta \eta_x$ median');
p4_4 = plot(t,slopeDiffMode,'LineWidth',2,'Color',co(3,:),'DisplayName', '$\Delta \eta_x$ mode');

ylabel(ax4_3,{'$\Delta \eta_x$ across troughs'},'Interpreter','latex')
ax4_3.YColor = [0,0,0];
set(ax4_3,'FontSize',18,'TickLabelInterpreter','latex','XLim',trange,'YLim',[0,0.4])
legend([p4_1,p4_2,p4_3,p4_4],'Interpreter','latex','Location', 'northwest')
linkaxes([ax1_1, ax4_3],'x')
%% Plot hovemollerish thing V2: don't load everything all at once
ti = 0;
dtPlot = 0.5;
tl = [ti, ti+dtPlot];
P_threshold = 7;
lambda_max = 4e-2; % asymmetry detection wavelength maximum.
slopePThreshold = 0.1;
slopeMin = [0.3,0.5,0.7];
solitonsThreshold = [0.01, 0.03]; % difference between distances to adjacent trough, minimum distance to an adjacent trough.

while true

    
    figure(2)
    
    x = (1:size(CompImg,2))*mpp;
    y = 1:size(CompImg,1);
    t = (y-1)/dydt;
    
    it = [0,0];
    [~, it(1)] = min(abs(t-tl(1)));
    [~, it(2)] = min(abs(t-tl(2)));
    
    
    hold off
    imagesc(CompImg(it(1):it(2),:),'XData',x,'YData',t(it(1):it(2)),[0,3])
    hold on
    
    %s = DataPath(end-23:end-1) + " " + frames(1) + " to " + frames(end);
    il = [0,0];
    il(1) = floor(tl(1)/spp)*2 + min(floor(mod(tl(1),spp)/dt_pair),1)+1;
    il(2) = floor(tl(2)/spp)*2 + min(floor(mod(tl(2),spp)/dt_pair),1)+1;
    
    s = sprintf('%s Frames %d to %d', [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)], frames(il(1)), frames(il(2)));
    title(s,'Interpreter','latex')
    xlabel('x (m)','Interpreter','latex')
    ylabel('t (s)','Interpreter','latex')
    
    % s = sprintf('%.4f (PN%.2f)',0:spp:t(end),1:t(end)/spp)
    set(gca,'DataAspectRatio',[1*mpp 1/dydt 1])
    set(gca, 'ytick', (round(t(it(1))/spp)*spp):spp:(round(t(it(2))/spp)*spp),'TickLabelInterpreter','latex','YTickLabel',compose("%.4f (PN %d)",((round(t(it(1))/spp)*spp):spp:(round(t(it(2))/spp)*spp))',(round(t(it(1))/spp):1:round(t(it(2))/spp))'));
    set(gca, 'xtick', 0:0.02:x(end));
    set(gca,'FontSize',24)
    colormap gray
    
    minPts = NaN(300,2);
    cmp = 1;
    for i = il(1):il(2)
        [TF, P] = islocalmin(-fYs(i,:));

        slopeMask = islocalmax(abs(diff(fYs(i,:))));
        steepPtsMask(1,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(1);
        steepPtsMask(2,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(2);
        steepPtsMask(3,:) = slopeMask & abs(diff(fYs(i,:))) >= slopeMin(3);
        
        minsMask = P > P_threshold;
        imins = find(minsMask);
        xmins = fXs(i,minsMask);
        ymins = fYs(i,minsMask);


        if ~isempty(ymins)
            slope = -diff(fYs(i,:));
            [slopeMins,slopeMinsP] = islocalmin(slope);
            [slopeMaxs,slopeMaxsP] = islocalmax(slope);

            if mod(i,2) == 1
                offset = round((i-1)/2*spp*dydt)+1-yLimits(1);
                offset2 = round((i-1)/2*spp*dydt)+1;
            else
                offset = round(((i-2)/2*spp+dt_pair)*dydt)+1-yLimits(1);
                offset2 = round(((i-2)/2*spp+dt_pair)*dydt)+1;
            end

            for j = 1:length(ymins)
                foundSlopeMin = false;
                foundSlopeMax = false;
                k = imins(j);
                troughSlopeMin = 0;
                troughSlopeMax = 0;
                while ~foundSlopeMin && k > 0 && imins(j)-k < lambda_max/2/mpp
                    if slopeMinsP(k) > slopePThreshold
                        foundSlopeMin = true;
                        troughSlopeMin = slope(k);
                        plot(fXs(i,k)*mpp,(fYs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                    end
                    k = k-1;
                end

                k = imins(j);
                while ~foundSlopeMax && k < length(slope) && k-imins(j) < lambda_max/2/mpp
                    if slopeMaxsP(k) > slopePThreshold
                        foundSlopeMax = true;
                        troughSlopeMax = slope(k);
                        plot(fXs(i,k)*mpp,(fYs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                    end
                    k = k+1;
                end

                slopeDiff = abs(troughSlopeMax + troughSlopeMin);
                slopeDiffText = sprintf('SlopeDiff %.2f',slopeDiff);
                text(xmins(j)*mpp,(ymins(j)+offset)/dydt,slopeDiffText)
                    
            end
        end
        
        solitonsMask = [false, abs(diff(diff(fXs(i,minsMask)))) > solitonsThreshold(1)/mpp, false];
        solitonsMask = solitonsMask | [false, diff(fXs(i,minsMask))>solitonsThreshold(2)/mpp] | [diff(fXs(i,minsMask))>solitonsThreshold(2)/mpp, false];
        if length(xmins) > 1
            solitonsMask(1) = solitonsMask(1) | xmins(1)>solitonsThreshold(2)/mpp;
            solitonsMask(end) = solitonsMask(end) | (fXs(end)-xmins(end) > solitonsThreshold(2)/mpp);

            solitonXmins = xmins(solitonsMask);
            solitonYmins = ymins(solitonsMask);
        elseif ~isempty(xmins)
            solitonsMask = true;
            solitonXmins = xmins(solitonsMask);
            solitonYmins = ymins(solitonsMask);
        end
    
        if mod(i,2) == 1
            offset = round((i-1)/2*spp*dydt)+1-yLimits(1);
            offset2 = round((i-1)/2*spp*dydt)+1;
        else
            offset = round(((i-2)/2*spp+dt_pair)*dydt)+1-yLimits(1);
            offset2 = round(((i-2)/2*spp+dt_pair)*dydt)+1;
        end

        plot(fXs(i,:)*mpp,(fYs(i,:)+offset)/dydt,'r')
        
        if ~isempty(xmins)
            plot(fXs(i,minsMask)*mpp,(fYs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
            plot(fXs(i,steepPtsMask(1,:))*mpp,(fYs(i,steepPtsMask(1,:))+offset)/dydt,'.g','MarkerSize',10)
            plot(fXs(i,steepPtsMask(2,:))*mpp,(fYs(i,steepPtsMask(2,:))+offset)/dydt,'.m','MarkerSize',20)
            plot(fXs(i,steepPtsMask(3,:))*mpp,(fYs(i,steepPtsMask(3,:))+offset)/dydt,'.y','MarkerSize',20)
            plot(solitonXmins*mpp,(solitonYmins+offset)/dydt,'.b','MarkerSize',10)
        end
    
        minPts(cmp:cmp+length(fXs(i,minsMask))-1,:) = [fXs(i,minsMask)',offset2*ones(1,length(fXs(i,minsMask)))'];
        cmp = cmp + length(fXs(i,minsMask));
        hold on
        % pause
    end
    minPts(any(isnan(minPts), 2), :) = [];
    
    ip = input('a for back, d for forward','s');
    nip = str2double(ip);
    scrollFrac = 0.5;
    if ip == 'a'
        tl(1) = max(0,tl(1)-dtPlot*scrollFrac);
        tl(2) = tl(1) + dtPlot;
    elseif ip == 'd'
        tl(2) = min(t(end),tl(2)+dtPlot*scrollFrac);
        tl(1) = tl(2) - dtPlot;
    elseif ~isnan(nip) && nip >= 0 && nip < t(end)-dtPlot
        tl(1) = nip;
        tl(2) = tl(1)+dtPlot;
    end
    tl
end
%% Plot hovemollerish thing
P_threshold = 10;


figure(2)

x = (1:size(CompImg,2))*mpp;
y = 1:size(CompImg,1);
t = y/dydt;

hold off
imagesc(CompImg,'XData',x,'YData',t,[0,3])
hold on

xlabel('x (m)','Interpreter','latex')
ylabel('t (s)','Interpreter','latex')

% s = sprintf('%.4f (PN%.2f)',0:spp:t(end),1:t(end)/spp)
set(gca,'DataAspectRatio',[1*mpp 1/dydt 1])
set(gca, 'ytick', 0:spp:t(end),'TickLabelInterpreter','latex','YTickLabel',compose("%.4f (PN %d)",(0:spp:t(end))',(0:1:round(t(end)/spp))'));
set(gca, 'xtick', 0:0.02:x(end));
set(gca,'FontSize',24)
colormap gray

ylim([29,30]);
tl = ylim;
%s = DataPath(end-23:end-1) + " " + frames(1) + " to " + frames(end);
il = [0,0];
il(1) = floor(tl(1)/spp)*2 + min(floor(mod(tl(1),spp)/dt_pair),1);
il(2) = floor(tl(2)/spp)*2 + min(floor(mod(tl(2),spp)/dt_pair),1);
s = sprintf('%s Frames %d to %d', [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)], frames(il(1)), frames(il(2)));
title(s,'Interpreter','latex')


minPts = NaN(300,2);
cmp = 1;
for i = 1:nF
    [TF, P] = islocalmin(-fYs(i,:));
    minsMask = P > P_threshold;

    if mod(i,2) == 1
        offset = round((i-1)/2*spp*dydt)+1-yLimits(1);
        offset2 = round((i-1)/2*spp*dydt)+1;
        plot(fXs(i,:)*mpp,(fYs(i,:)+offset)/dydt,'r')
        plot(fXs(i,minsMask)*mpp,(fYs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
    else
        offset = round(((i-2)/2*spp+dt_pair)*dydt)+1-yLimits(1);
        offset2 = round(((i-2)/2*spp+dt_pair)*dydt)+1;
        plot(fXs(i,:)*mpp,(fYs(i,:)+offset)/dydt,'r')
        plot(fXs(i,minsMask)*mpp,(fYs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
    end
    
    minPts(cmp:cmp+length(fXs(i,minsMask))-1,:) = [fXs(i,minsMask)',offset2*ones(1,length(fXs(i,minsMask)))'];
    cmp = cmp + length(fXs(i,minsMask));
    hold on
    % pause
end
minPts(any(isnan(minPts), 2), :) = [];


% plot(minPts(:,1),minPts(:,2),'r*') 
%%
troughImg = zeros(size(CompImg));
idxs = sub2ind(size(troughImg), minPts(:,2),minPts(:,1));
troughImg(idxs) = 1;
troughImg = imgaussfilt(troughImg,10);
figure(5)

troughImg = (troughImg > 0.00001);
imagesc(troughImg)

figure(3)
[H,theta,rho]= hough(troughImg,'Theta',-45:0.01:-10);
% Hblur = imgaussfilt(H,5);
imshow(imadjust(rescale(H)),'XData',theta,'YData',rho,...
      'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);

P_threshold = 0.15;
[TF1, P1] = islocalmax(Hblur,1);
[TF2, P2] = islocalmax(Hblur,2);
P = P1.*P2;
figure(4)
imagesc(P,'XData',theta,'YData',rho)
[iPr,iPt] = find(P > P_threshold);
thetas = theta(iPt);
rhos = rho(iPr);

figure(1)
hold on
ms = NaN(1,length(thetas));
for i = 1:length(thetas)
    x = get(gca,'XLim');
    ms(i) = -cos(theta(i))/sin(theta(i));
    y = (rhos(i) - x* cos(theta(i)))/ sin(theta(i));
    plot(x,y,'r');
end
% plot(XPIVSurfW1_Surface, PIVSurfW1_Surface, 'r')
%%
% plotLine([1734,1366], [1978,1571],1)
plotLine([2091,1972],[2335,2169],1)
plotLine([242,1963],[477,2163],1)
plotLine([1621,923],[1388,725],1)
plotLine([327,1341],[574,1544],1)

%%
function plotLine(A,B,figNum)
    figure(figNum)
    xlim = get(gca,'XLim');
    m = (B(2)-A(2))/(B(1)-A(1));
    b = B(2) - m*B(1);
    y1 = m*xlim(1) + b;
    y2 = m*xlim(2) + b;
    hold on
    plot([xlim(1) xlim(2)],[y1 y2],'-r','LineWidth',2)
    plot([A(1) B(1)],[A(2) B(2)],'*r','MarkerSize',20)
    hold off
end
%%
function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface] = FindWaterSurface(PIVSurfW_CamAngle)
    X = 31:size(PIVSurfW_CamAngle,2)-40;
    [imSurf] = Copy_of_FindSurface(PIVSurfW_CamAngle(1:2800,X), 5, 5);
    PIVSurf_Surface_Raw = imSurf.surface;
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
    PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
    [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1d3);
    if length(SP.coefs)>2
        PIVSurf_Surface_W = SP.coefs(2:end-1);
    else
        CC = polyfit([X(1),X(end)],[SP.coefs(1) SP.coefs(2)],1);
        PIVSurf_Surface_W = polyval(CC,[X(1):X(end)]);
    end
    Usurf = X;
    Vsurf = PIVSurf_Surface_W;
    
    %%% Check if bad frame
    BadFramePIVSurfW = 0;
    if imSurf.badFrameBool == 1
        BadFramePIVSurfW = 1;
    end
    
    XPIVSurfW_Surface = Usurf;
    PIVSurfW_Surface = Vsurf;
end