clear
clc

%%

LoadPath = '/media/surflab/Working24/ExpAW/ExpAW5_acc0.22_W5V/ExpAW5_acc0.22_W5V_Run2/RESULTS_andy/Air/';

PIVWaterDir = dir([LoadPath 'PIV_Velocities_raw/' '*.mat']); %Same for water

PIV1Dir_temp = PIVWaterDir;

mpp = 6.493178e-5; %Accurate only for the main dataset, not the pilot

% -ExpAW1Run2 around 1340 but not as prominant
% -ExpAW5 Run2 super steep parasitic? capillaries. Perhaps periodic crapper
% capillaries. solitons? starting around 824. Steep soliton around 861.
% Some good asymmetry
%% Save frames for a video


clear f F


figure('units','pixels','Position',[0,0,1000,1000])

f = 1;
tic

idxs = 1:length(PIV1Dir_temp);
% F = struct('cdata',cell(length(idxs),1),'colormap',cell(length(idxs),1));
parfor i = 1:length(idxs)

    idx = idxs(i);
    fname = [PIV1Dir_temp(idx).folder '/' PIV1Dir_temp(idx).name];
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