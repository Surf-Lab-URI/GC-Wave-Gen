clear
clc

%% Load Raw Velocities

LoadPath = '/media/surflab/Working24/ExpAW/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/RESULTS_andy/Air/';

PIVAirDir_RAW = dir([LoadPath 'PIV_Velocities_raw/' '*.mat']); %Same for water

%% Save Raw Velocity frames for a video


clear f F


figure('units','pixels','Position',[0,0,1000,1000])

f = 1;
tic

idxs = 1:length(PIVAirDir_RAW);
% F = struct('cdata',cell(length(idxs),1),'colormap',cell(length(idxs),1));
parfor i = 1:length(idxs)

    idx = idxs(i);
    fname = [PIVAirDir_RAW(idx).folder '/' PIVAirDir_RAW(idx).name];
    CompVelAir = load(fname);
    CST = CompVelAir.CST;

    hold off
    imagesc(CompVelAir.delta_x*CST.DX/CST.DT,[0,6])
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar;
    c.Label.String = "Horizontal Velocity (m/s)";
    colormap gray
    axis off
    set(gca,'FontSize',24)
    drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb)
        delete(sbt)
    end
    sb = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    drawnow

    fname = ['videoframes/' LoadPath(end-40:end-18) '_' num2str(CompVelAir.PairNum) '.jpg']; % full name of image
    % print('-djpeg','-r600',fname)     % save image with '-r200' resolution
    % saveas(gcf,fname,'tiffn')
    saveas(gcf,fname)
end
toc



%% Load smoothed velocities
clear

LoadPath = '/media/surflab/Working24/ExpAW/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/RESULTS_andy/Air/';

PIVAirDir_Cart = dir([LoadPath 'CALCULATED_FIELDS/Cartesian Fields/Velocity/' '*.mat']); %Same for water

SaveUImgDir_Cart = [LoadPath 'CALCULATED_FIELDS/Cartesian Fields/Velocity/U/'];
SaveWImgDir_Cart = [LoadPath 'CALCULATED_FIELDS/Cartesian Fields/Velocity/W/'];
SaveSpeedImgDir_Cart = [LoadPath 'CALCULATED_FIELDS/Cartesian Fields/Velocity/Speed/'];


if ~exist(SaveUImgDir_Cart, 'dir')
    mkdir(SaveUImgDir_Cart);
end

if ~exist(SaveWImgDir_Cart, 'dir')
    mkdir(SaveWImgDir_Cart);
end

if ~exist(SaveSpeedImgDir_Cart, 'dir')
    mkdir(SaveSpeedImgDir_Cart);
end

%% Save Smoothed Velocity frames for a video
clear f

velfig = figure('units','pixels','Position',[0,0,1200,2400]);
drawnow
f = 1;
tic

qdx = 20;

idxs = 1:length(PIVAirDir_Cart);

parfor i = 1:length(idxs)

    idx = idxs(i);
    fname = [PIVAirDir_Cart(idx).folder '/' PIVAirDir_Cart(idx).name];
    Cartesian_Air = load(fname);
    CST = Cartesian_Air.CST;

    u = Cartesian_Air.u.*Cartesian_Air.Mask*CST.DX/CST.DT;
    w = Cartesian_Air.w.*Cartesian_Air.Mask*CST.DX/CST.DT; %Positive is up for w and z. Postive is down for v and y.
    speed = (u.^2 + w.^2).^0.5;
    mask = Cartesian_Air.Mask;
    
    %% Plot U
    
    hold off
    T = tiledlayout(4,1,'TileSpacing','tight');
    ax1 = nexttile(2,[3 1]);
    qdx = 20;
    cl = [0,6];
    imagesc(u,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar; 
    c.Label.String = "Horizontal Velocity (m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[1,630])
    
    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.03)*[1 1];
    try
        delete(sb1)
        delete(sbt1)
    end
    sb1 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt1 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    % drawnow



    ax2 = nexttile(1);
    hold off
    qdx = 10;
    cl = [0,4];
    imagesc(u,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar;
    c.Label.String = "(m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[430,580])
    linkaxes([ax1 ax2], 'x')
    % ax2.YLim = [430,550];

    

    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb2)
        delete(sbt2)
    end
    sb2 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt2 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    drawnow

    fname = [SaveUImgDir_Cart LoadPath(end-40:end-18) '_' num2str(Cartesian_Air.PairNum) '_U_Cart.svg']; % full name of image
    saveas(gcf,fname)

    %% Plot W
    
    hold off
    T = tiledlayout(4,1,'TileSpacing','tight');
    ax1 = nexttile(2,[3 1]);
    qdx = 20;
    cl = [-1,1];
    imagesc(w,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar; 
    c.Label.String = "Vertical Velocity (m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[1,630])
    
    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.03)*[1 1];
    try
        delete(sb1)
        delete(sbt1)
    end
    sb1 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt1 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    % drawnow



    ax2 = nexttile(1);
    hold off
    qdx = 10;
    cl = [-1,1];
    imagesc(u,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar;
    c.Label.String = "(m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[430,580])
    linkaxes([ax1 ax2], 'x')
    % ax2.YLim = [430,550];

    

    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb2)
        delete(sbt2)
    end
    sb2 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt2 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    drawnow

    fname = [SaveWImgDir_Cart LoadPath(end-40:end-18) '_' num2str(Cartesian_Air.PairNum) '_W_Cart.jpg']; % full name of image
    saveas(gcf,fname)


    %% Plot Speed
    
    hold off
    T = tiledlayout(4,1,'TileSpacing','tight');
    ax1 = nexttile(2,[3 1]);
    qdx = 20;
    cl = [0,6];
    imagesc(speed,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar; 
    c.Label.String = "Speed (m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[1,630])
    
    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.03)*[1 1];
    try
        delete(sb1)
        delete(sbt1)
    end
    sb1 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt1 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    % drawnow



    ax2 = nexttile(1);
    hold off
    qdx = 10;
    cl = [0,4];
    imagesc(u,cl)
    hold on
    set(gca,'DataAspectRatio',[1 1 1])
    c = colorbar;
    c.Label.String = "(m/s)";
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    colormap gray
    axis off
    set(gca,'FontSize',24,'YLim',[430,580])
    linkaxes([ax1 ax2], 'x')
    % ax2.YLim = [430,550];

    

    %Note: plotting quiver on top of imagesc inverts vertical axis.
    quiver(1:qdx:size(mask,2), 1:qdx:size(mask,1),u(1:qdx:end,1:qdx:end), -w(1:qdx:end,1:qdx:end),'r')

    % drawnow

    xl = xlim;
    yl = ylim;
    
    lsbm = 1e-2; % length of scale bar in meters
    lsb = lsbm/CST.DX/CST.GS;
    xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
    ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
    try
        delete(sb2)
        delete(sbt2)
    end
    sb2 = plot(xsb,ysb,'-r', 'LineWidth',5);
    
    
    sbl = sprintf('%d cm',lsbm*100);
    sbt2 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');
    drawnow

    fname = [SaveSpeedImgDir_Cart LoadPath(end-40:end-18) '_' num2str(Cartesian_Air.PairNum) '_Speed_Cart.jpg']; % full name of image
    saveas(gcf,fname)

end
toc
close(velfig)

%% Calculate u'w' over time
clear


LoadPath = '/media/surflab/Working24/ExpAW/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/RESULTS_andy/Air/';

PIVAirDir_Cart = dir([LoadPath 'CALCULATED_FIELDS/Cartesian Fields/Velocity/' '*.mat']); %Same for water

idx = 1;
fname = [PIVAirDir_Cart(idx).folder '/' PIVAirDir_Cart(idx).name];
Cartesian_Air = load(fname);
CST = Cartesian_Air.CST;

up_b = zeros(length(PIVAirDir_Cart),size(Cartesian_Air.Mask,1));
wp_b = zeros(length(PIVAirDir_Cart),size(Cartesian_Air.Mask,1));
u_b = zeros(length(PIVAirDir_Cart),1);
w_b = zeros(length(PIVAirDir_Cart),1);
upwp_b = zeros(length(PIVAirDir_Cart),size(Cartesian_Air.Mask,1));
dupwpdz_b = zeros(length(PIVAirDir_Cart),size(Cartesian_Air.Mask,1)-1);
PairNums = zeros(length(PIVAirDir_Cart),1);

idxs = 1:length(PIVAirDir_Cart);

tic
parfor i = 1:length(idxs)

    idx = idxs(i);
    fname = [PIVAirDir_Cart(idx).folder '/' PIVAirDir_Cart(idx).name];
    Cartesian_Air = load(fname);
    CST = Cartesian_Air.CST;

    u = Cartesian_Air.u.*Cartesian_Air.Mask*CST.DX/CST.DT;
    w = Cartesian_Air.w.*Cartesian_Air.Mask*CST.DX/CST.DT; %Positive is up for w and z. Postive is down for v and y.
    speed = (u.^2 + w.^2).^0.5;
    mask = Cartesian_Air.Mask;
    
    u_b = mean(u,2,'omitnan');
    w_b = mean(w,2,'omitnan');
    
    up = u - u_b;
    wp = w - w_b;

    up_b(i,:) = mean(up,2,'omitnan');
    wp_b(i,:) = mean(wp,2,'omitnan');
    
    upwp_b(i,:) = mean(up.*wp,2,'omitnan');
    
    dupwpdz_b(i,:) = mean(diff(up.*wp,1,1),2,'omitnan');

    PairNums(i) = str2double(Cartesian_Air.PairNum);
    disp(['Pair ' Cartesian_Air.PairNum ' Finished!'])

end
toc

pps = 14.5; % pairs per second
spp = 1/pps; % seconds per pair
t = PairNums*spp;

figure
plot(PairNums, movmean(upwp_b(:,500),10))
%% Plot u'w'
figure('units','pixels','Position',[0,0,1000,1000])
plot(PairNums, movmean(upwp_b(:,500),30),'LineWidth',4)
hold on
xlabel('Pair Number','Interpreter','latex')
ylabel("$\overline{u'v'}\ \mathrm{(m^2/s^2)}$",'Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')
title('0.6 cm above mean water level, 30 pair moving mean','Interpreter','latex')

figure('units','pixels','Position',[0,0,1000,1000])
plot(t, movmean(upwp_b(:,500),30),'LineWidth',4)
hold on
xlabel('Time (s)','Interpreter','latex')
ylabel("$\overline{u'v'}\ \mathrm{(m^2/s^2)}$",'Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')
title('0.6 cm above mean water level, 30 pair moving mean','Interpreter','latex')

%%
figure('units','pixels','Position',[0,0,1000,1000])
imagesc(movmean(upwp_b,10,1)','XData',t,[-0.05,0.02])
hold on
colormap gray
xlim([0,50])
ylim([0,600])
set(gca,'ytick',[])
xlabel('Time (s)','Interpreter','latex')
set(gca,'FontSize',30)
set(gca,'TickLabelInterpreter','latex')
c = colorbar;
c.Label.String = "$\overline{u'v'}\ \mathrm{(m^2/s^2)}$";
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
drawnow

xl = xlim;
yl = ylim;

lsbm = 1e-2; % length of scale bar in meters
lsb = lsbm/CST.DX/CST.GS;
xsb = ((xl(2)-xl(1))*0.03)*[1 1];
ysb = [yl(2)-(yl(2)-yl(1))*0.12, yl(2)-(yl(2)-yl(1))*0.12+lsb];
try
    delete(sb1)
    delete(sbt1)
end
sb1 = plot(xsb,ysb,'-r', 'LineWidth',5);


sbl = sprintf('%d cm',lsbm*100);
sbt1 = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'Color','red','FontSize',16,'Interpreter','latex');



%% Compute surface height and slope variance
% clear

pps = 14.5; % pairs per second

spp = 1/pps; % seconds per pair

dt_pair = 22.222e-3; % time between pictures in a given pair

LoadPath = '/media/surflab/Working24/ExpAW/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/RESULTS_andy/';

WaterSurfDir = dir([LoadPath 'Water/Surfaces/' '*.mat']);
nPairs = length(WaterSurfDir);
nF = nPairs*2;

%Load first pair to get dimensions
idx = 1;
fname = [WaterSurfDir(idx).folder '/' WaterSurfDir(idx).name];
SavedSurfsWater = load(fname);
CST = SavedSurfsWater.CST;
PixRes_Water1 = SavedSurfsWater.PixRes_Water1;
PixRes_Water2 = SavedSurfsWater.PixRes_Water2;

fXs = zeros(nF, length(PixRes_Water1.PIVW_PIVSurfW1_Surface));
fYs = zeros(nF, length(PixRes_Water1.PIVW_PIVSurfW1_Surface));
pairIndices = zeros(nF,1);
pairnums = zeros(nF,1);
ImageNums = zeros(nF,1);
t = zeros(nF,1);

tic
for i = 1:length(WaterSurfDir)
    fname = [WaterSurfDir(i).folder '/' WaterSurfDir(i).name];
    SavedSurfsWater = load(fname);
    pair_index = SavedSurfsWater.PixRes_Water1.pair_index;
    

    PixRes_Water1 = SavedSurfsWater.PixRes_Water1;
    PixRes_Water2 = SavedSurfsWater.PixRes_Water2;
    
    fXs(i*2-1,:) = PixRes_Water1.XPIVW_PIVSurfW1_Surface;
    fYs(i*2-1,:) = PixRes_Water1.PIVW_PIVSurfW1_Surface;
    pairIndices(i*2-1) = PixRes_Water1.pair_index;
    pairnums(i*2-1) = str2num(PixRes_Water1.PairNum);
    ImageNums(i*2-1) = str2num(PixRes_Water1.ImageNum_1);
    t(i*2-1) = pairnums(i*2-1)*spp;

    fXs(i*2,:) = PixRes_Water2.XPIVW_PIVSurfW2_Surface;
    fYs(i*2,:) = PixRes_Water2.PIVW_PIVSurfW2_Surface;
    pairIndices(i*2) = PixRes_Water2.pair_index;
    pairnums(i*2) = str2num(PixRes_Water1.PairNum);
    ImageNums(i*2) = str2num(PixRes_Water1.ImageNum_1);
    t(i*2) = pairnums(i*2)*spp + dt_pair;
end
toc

% figure(2)
% for i = 1:nF
%     plot(fXs(i,:), -fYs(i,:))
%     set(gca,'DataAspectRatio',[1 1 1])
%     ylim([-800,-600])
%     pause(0.1)
% end

mpp = CST.DX_W;
eta = (fYs-mean(fYs(1:20,:),'all'))*mpp;
x_eta = fXs*mpp;
eta_var = sum((eta-mean(eta,2)).^2,2)/(size(eta,2)-1);

eta_x = diff(eta/mpp,1,2);

eta_x_var = sum((eta_x-mean(eta_x,2)).^2,2)/(size(eta_x,2)-1);

figure(3)
ax1 = subplot(3,1,1);
plot(t,eta_var,'LineWidth',4)
hold on
xlabel('Time (s)','Interpreter','latex')
ylabel('$\mathrm{Var}[\eta]\ \mathrm{(m^2)}$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')


ax2 = subplot(3,1,2);
plot(t,eta_x_var,'LineWidth',4)
hold on
xlabel('Time (s)','Interpreter','latex')
ylabel('$\mathrm{Var}[\eta_x]$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')
xlim([30,50])
ylim([0,0.1])

%Calculate Surface Area
A_surf = zeros(1,length(t));
for n = 1:length(t)
    for i = 1:(size(eta,2)-1)
        A_surf(n) = A_surf(n) + ((eta(n,i+1) - eta(n,i))^2 + (x_eta(n,i+1) - x_eta(n,i))^2).^0.5;
    end
end


ax3 = subplot(3,1,3);
l_interface = x_eta(1,end) - x_eta(1,1);
plot(t,(A_surf - l_interface)/l_interface*100,'LineWidth',4)
% plot(t,A_surf,'LineWidth',4)
hold on
xlabel('Time (s)','Interpreter','latex')
ylabel('$\mathrm{Surface\ Area\ Increase\ (\%)}$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')
linkaxes([ax1,ax2,ax3],'x')
xlim([30,50])
ylim([0,4])

figure(5)
plot(ImageNums,eta_x_var)
hold on
xlabel('Frame Number')
ylabel('Spacial variance of surface slope at each time (m/m)^2')
set(gca,'FontSize',24)

%% Generate Video File from frames save in previous section
vw = VideoWriter('ExpAW5_acc0.22_W5V_Run2_750to1350_ED.avi', 'Uncompressed AVI');
vw.FrameRate = 2;
open(vw);
vw

for i = idxs(1):2:idxs(end-1)
    fname = ['videoframes/' DataPath(end-23:end-1) '_' num2str(i)]; % full name of image
    I = imread([fname '.jpg']); % read saved image
    frame = im2frame(I(4:1190,10:8685,:)); % convert image to frame
    writeVideo(vw,frame)
end

close(vw);