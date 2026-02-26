clear
clc
close all

LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';

resultsfname = [LONG, 'Results_Surflab/results.mat'];
load(resultsfname,"exps")

%% Plot surface variance for one experiment
ii = 4;
exp_name = exps{ii}.exp_name
Surfs = exps{ii}.Surfs;
eta = Surfs.eta;
surfs = Surfs.surfs;
x_eta = (0:(size(eta,2)-1))*Surfs.dx;
dx = Surfs.dx;
dt_pair = Surfs.dt_pair;
t_eta = Surfs.t;
nP = exps{ii}.number_of_pair;
mean(eta,'all')
mean(eta(1:20),"all")
mean(surfs*dx,'all')
mean(surfs(1:20,:)*dx,'all')
eta_var = sum((eta-mean(eta,2)).^2,2)/(size(eta,2)-1);

trange = [10,40];

figure(1)
plot(t_eta, eta_var)
set(gca, 'YLim',[0,2.5e-7],'XLim',trange)
%% Plot IR cam temperature and surface elevation together
IRDir = dir([LONG exp_name '/IRMat/']);
IRDir = IRDir(3:end-1);
itemp = 1;
clear tempss eta_PIV_xs
for tImg = 17.8:(1/7.2):21
[~,surfIdx] = min(abs(t_eta-tImg));
if mod(surfIdx,2) == 0
    surfIdx = surfIdx-1;
end
surfIdx;
t_eta(surfIdx)

pairNum = Surfs.pairNum(surfIdx);

IRNum = pairNum*6;

if str2double(IRDir(IRNum+1).name(end-7:end-4))~=IRNum
    disp("Warning: IR images missing. Number of file name doesn't match position in directory. Results may be wrong")
end

IRPath = IRDir(IRNum+1);

load([IRPath.folder,'/',IRPath.name],'IR');

figure(1)
imagesc(IR.img,[19.5,20.5])
colorbar
hold on
plot([290,290],[202,511],'-r','LineWidth',2)

ir = IR.img(202:511,290)';
x_ir = IR.DX*(0:length(ir)-1)-0.002;

eta_PIV = -CropSurfToPIVDims(eta(surfIdx,:),false);
x_eta_PIV = (0:length(eta_PIV)-1)*Surfs.dx;

figure(2)
ax1 = subplot(2,1,1);
hold off
plot(x_eta_PIV,eta_PIV,'LineWidth',3,'DisplayName','$\eta$')
hold on
daspect([1,1,1])
ylim([-0.005,0.005])
ax2 = subplot(2,1,2);
hold off
plot(x_ir,ir,'LineWidth',2,'DisplayName','Temperature')
hold on
xlim([x_ir(1),x_ir(end)])
linkaxes([ax1,ax2],'x')

eta_PIV_x = diff(eta_PIV)/Surfs.dx;
temps = nan(1,length(eta_PIV_x));
for i = 1:length(temps)
    [~,idx_ir] = min(abs(x_eta_PIV(i)-x_ir));
    temps(i) = ir(idx_ir);
end
eta_PIV_xs(itemp:(itemp+length(eta_PIV_x)-1)) = eta_PIV_x;
tempss(itemp:(itemp+length(eta_PIV_x)-1)) = temps;
itemp = itemp + length(eta_PIV_x);
figure(3)
plot(temps,abs(eta_PIV_x),'.r','MarkerSize',5)
hold on
set(gca,'FontSize',24,'TickLabelInterpreter','latex')
xlabel('T ($\circ$ C)','Interpreter','latex')
ylabel('$\eta_x$', 'Interpreter','latex')
corrcoef(tempss, abs(eta_PIV_xs))

pause
end



%% Plot hovemollerish thing V2: don't load everything all at once
% Assemble composite image for hovmoller plot
tic

nF = nP*2;

ymean = 1790;
yLimits = [-1,1]*30 + ymean;
dy = yLimits(2)-yLimits(1);

dydt = dy/dt_pair;

spp = Surfs.spp;

DeltaT = nP*spp;

CompImg = NaN(ceil(DeltaT*dydt),3308);
for idx = 1:nF
    imgPath = Surfs.paths(idx);
    imgPath
    
    load(imgPath,"imgPivsurf");

    [~,surfImg,~] = SurfImgToPIVDims(imgPivsurf);
     
    if mod(idx,2) == 1
        CompImg(round((idx-1)/2*spp*dydt)+1:round((idx-1)/2*spp*dydt)+dy+1,1:size(surfImg,2)) = surfImg(yLimits(1):yLimits(2),:)/mean(surfImg,'all');
    else
        CompImg(round(((idx-2)/2*spp+dt_pair)*dydt)+1:round(((idx-2)/2*spp+dt_pair)*dydt)+dy+1,1:size(surfImg,2)) = surfImg(yLimits(1):yLimits(2),:)/mean(surfImg,'all');
    end
end
toc
%% Plot interactive Hovmollerish thing

ti = 0;
dtPlot = 0.5;
tl = [ti, ti+dtPlot];
P_threshold = 7;
lambda_max = 4e-2; % asymmetry detection wavelength maximum.
slopePThreshold = 0.1;
slopeMin = [0.3,0.5,0.7];
solitonsThreshold = [0.01, 0.03]; % difference between distances to adjacent trough, minimum distance to an adjacent trough.
x_pix = 1:size(surfs,2);

while true
    
    figure(2)
    x = (1:size(CompImg,2))*dx;
    y = 1:size(CompImg,1);
    t = (y-1)/dydt;
    
    it = [0,0];
    [~, it(1)] = min(abs(t-tl(1)));
    [~, it(2)] = min(abs(t-tl(2)));
    
    
    hold off
    imagesc(CompImg(it(1):it(2),:),'XData',x,'YData',t(it(1):it(2)),[0.75 1.6])
    hold on
    
    il = [0,0];
    il(1) = floor(tl(1)/spp)*2 + min(floor(mod(tl(1),spp)/dt_pair),1)+1;
    il(2) = floor(tl(2)/spp)*2 + min(floor(mod(tl(2),spp)/dt_pair),1)+1;
    
    s = sprintf('%s Pairs %d to %d', exp_name, Surfs.pairNum(il(1)), Surfs.pairNum(il(2)));
    title(s,'Interpreter','latex')
    xlabel('x (m)','Interpreter','latex')
    ylabel('t (s)','Interpreter','latex')
    
    % s = sprintf('%.4f (PN%.2f)',0:spp:t(end),1:t(end)/spp)
    set(gca,'DataAspectRatio',[1*dx 1/dydt 1])
    set(gca, 'ytick', (round(t(it(1))/spp)*spp):spp:(round(t(it(2))/spp)*spp),'TickLabelInterpreter','latex','YTickLabel',compose("%.4f (PN %d)",((round(t(it(1))/spp)*spp):spp:(round(t(it(2))/spp)*spp))',(round(t(it(1))/spp):1:round(t(it(2))/spp))'));
    set(gca, 'xtick', 0:0.02:x(end));
    set(gca,'FontSize',24)
    colormap gray
    
    minPts = NaN(300,2);
    cmp = 1;
    for i = il(1):il(2)
        [TF, P] = islocalmin(-surfs(i,:));

        slopeMask = islocalmax(abs(diff(surfs(i,:))));
        steepPtsMask(1,:) = slopeMask & abs(diff(surfs(i,:))) >= slopeMin(1);
        steepPtsMask(2,:) = slopeMask & abs(diff(surfs(i,:))) >= slopeMin(2);
        steepPtsMask(3,:) = slopeMask & abs(diff(surfs(i,:))) >= slopeMin(3);
        
        minsMask = P > P_threshold;
        imins = find(minsMask);
        xmins = x_pix(minsMask);
        ymins = surfs(i,minsMask);


        if ~isempty(ymins)
            slope = -diff(surfs(i,:));
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
                while ~foundSlopeMin && k > 0 && imins(j)-k < lambda_max/2/dx
                    if slopeMinsP(k) > slopePThreshold
                        foundSlopeMin = true;
                        troughSlopeMin = slope(k);
                        plot(x_pix(k)*dx,(surfs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                    end
                    k = k-1;
                end

                k = imins(j);
                while ~foundSlopeMax && k < length(slope) && k-imins(j) < lambda_max/2/dx
                    if slopeMaxsP(k) > slopePThreshold
                        foundSlopeMax = true;
                        troughSlopeMax = slope(k);
                        plot(x_pix(k)*dx,(surfs(i,k)+offset)/dydt,'.k','MarkerSize',20)
                    end
                    k = k+1;
                end

                slopeDiff = abs(troughSlopeMax + troughSlopeMin);
                slopeDiffText = sprintf('SlopeDiff %.2f',slopeDiff);
                text(xmins(j)*dx,(ymins(j)+offset)/dydt,slopeDiffText)
                    
            end
        end
        
        solitonsMask = [false, abs(diff(diff(x_pix(minsMask)))) > solitonsThreshold(1)/dx, false];
        solitonsMask = solitonsMask | [false, diff(x_pix(minsMask))>solitonsThreshold(2)/dx] | [diff(x_pix(minsMask))>solitonsThreshold(2)/dx, false];
        if length(xmins) > 1
            solitonsMask(1) = solitonsMask(1) | xmins(1)>solitonsThreshold(2)/dx;
            solitonsMask(end) = solitonsMask(end) | (x_pix(end)-xmins(end) > solitonsThreshold(2)/dx);

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

        plot(x_pix(:)*dx,(surfs(i,:)+offset)/dydt,'r')
        
        if ~isempty(xmins)
            plot(x_pix(minsMask)*dx,(surfs(i,minsMask)+offset)/dydt,'.r','MarkerSize',10)
            plot(x_pix(steepPtsMask(1,:))*dx,(surfs(i,steepPtsMask(1,:))+offset)/dydt,'.g','MarkerSize',10)
            plot(x_pix(steepPtsMask(2,:))*dx,(surfs(i,steepPtsMask(2,:))+offset)/dydt,'.m','MarkerSize',20)
            plot(x_pix(steepPtsMask(3,:))*dx,(surfs(i,steepPtsMask(3,:))+offset)/dydt,'.y','MarkerSize',20)
            plot(solitonXmins*dx,(solitonYmins+offset)/dydt,'.b','MarkerSize',10)
        end
    
        minPts(cmp:cmp+length(x_pix(minsMask))-1,:) = [x_pix(minsMask)',offset2*ones(1,length(x_pix(minsMask)))'];
        cmp = cmp + length(x_pix(minsMask));
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
