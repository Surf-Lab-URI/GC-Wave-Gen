%% Step through raw frames
clear
clc
close all

% Define Path
DataPath = '/media/surflab/Working24/ExpAW/ExpAW4_acc0.16_W5V/ExpAW4_acc0.16_W5V_Run2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_Andy/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water
%%
frames = 0:1:1500;
nF = length(frames);
fXs = zeros(nF, 3639);
fYs = zeros(nF, 3639);

mpp = 6.493178e-5; %for main dataset only

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

%%
framesCtr = 1;
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
end
toc
save('temp.mat','-v7.3')
clear
load('temp.mat')
%%
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
%% Plot variance of eta, eta_x. Plot surface area
t = zeros(1,nF);
t(1) = floor(frames(1)/2)*spp + mod(frames(1),2)*dt_pair;
for i = 2:nF
    if mod(frames(i),2)==1
        t(i) = t(i-1)+dt_pair;
    else
        t(i) = t(i-1)+spp-dt_pair;
    end
end

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
ylim([0,2.5e-7])
s = [DataPath(end-23:end-18), '\_', DataPath(end-16:end-10), '\_', DataPath(end-8:end-6), '\_', DataPath(end-4:end-1)];
title(s,'Interpreter','latex')


ax2 = subplot(3,1,2);
plot(t,eta_x_var,'LineWidth',4)
hold on
xlabel('Time (s)','Interpreter','latex')
ylabel('$\mathrm{Var}[\eta_x]$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')
xlim([26,38])
ylim([0,0.05])

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
xlim([24,38])
ylim([0,4])


%% Fit growth rate to initial wavelet growth stage
%For ExpAW5R2
% tfit_i = 27;
% tfit_f = 29.58;

%For ExpAW4R2
tfit_i = 32.5;
tfit_f = 35.5;

ir = [0,0];
[~,i] = min(abs(t-tfit_i));
ir(1) = i;
[~,i] = min(abs(t-tfit_f)); 
ir(2) = i;
ir
eta_var_fit = fit(t(ir(1):ir(2))',eta_var(ir(1):ir(2)),'exp1')
% figure
% plot(eta_var_fit,t(ir(1):ir(2))',eta_var(ir(1):ir(2)))

figure(3)
hold on
subplot(3,1,1)
s = sprintf('Points %.2f s to %.2f s',tfit_i,tfit_f);
p1 = plot(t(ir(1):ir(2))',eta_var(ir(1):ir(2)),'.k', 'MarkerSize',15,'DisplayName',s)

fit_line_t = (t(ir(1)):0.01:t(ir(2)))';
fit_line = eta_var_fit.a.*exp(fit_line_t.*eta_var_fit.b);
s = sprintf('Fit %.2f s to %.2f s, Growth Rate %.2f 1/s',tfit_i,tfit_f,eta_var_fit.b);
p2 = plot(fit_line_t, fit_line, '-r', 'LineWidth',3,'DisplayName',s)
legend([p1,p2],'Interpreter','latex')

%% Plot Wavenumber Spectrum
nX = size(eta,2);
ppm = 1/mpp;

kspec_mag = abs(fft(eta-mean(eta,2),nX,2));

%Calculate contribution to total variance by each wavenumber, normalized by wavenumber bandwidth
%Sum of kspec_mag_n*(k step size) = kspec_mag_n*(2*pi/nX/mpp) is equal to
%the tocal variance
kspec_vareta_n = kspec_mag.^2/(nX-1)*mpp/(2*pi); 

figure(13)
hold off
imagesc(kspec_mag.^2/(nX-1)*mpp/(2*pi),'XData',2*pi*(ppm/nX*(0:nX-1)),'YData',t,[0,1e-9])
hold on
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlim([0,500])
ylim([32,36])
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


set(gca,'DataAspectRatio',[1*mpp 1/dydt 1])
set(gca, 'ytick', 0:spp:t(end),'TickLabelInterpreter','latex');
set(gca, 'xtick', 0:0.02:x(end));
set(gca,'FontSize',24)
colormap gray

ylim([34,35]);
yl = ylim;
%s = DataPath(end-23:end-1) + " " + frames(1) + " to " + frames(end);
il = [0,0];
il(1) = floor(yl(1)/spp)*2 + min(floor(mod(yl(1),spp)/dt_pair),1);
il(2) = floor(yl(2)/spp)*2 + min(floor(mod(yl(2),spp)/dt_pair),1);
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