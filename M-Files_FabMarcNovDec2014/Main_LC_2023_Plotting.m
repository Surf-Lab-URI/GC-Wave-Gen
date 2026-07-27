clear
% clc
LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1%:length(DIRS)

exp_name=DIRS(ii).name

num_of_digits = 3;
load_path = [LONG exp_name];
files=dir([load_path '/PIVRaw/PIV/*.mat']);
number_of_pair=length(files)/2;

%%%%%%%% Quickly flip through frames if you want (set previewing to true)

image_pair_number = 100;
previewing = false;
while image_pair_number < number_of_pair-1 && previewing
    load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
    IM_a = imgPiv;
    load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
    IM_b = imgPiv;

    %PIV Surf
    load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
    imgPivsurfa = imgPivsurf;
    load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
    imgPivsurfb = imgPivsurf;
    
    disp(['image_pair_number = ', num2str(image_pair_number)])
    figure(1)
    imagesc(imgPivsurfa, [0,300])
    colormap gray
    figure(2)
    imagesc(imgPivsurfb, [0,300])
    colormap gray
    
    ip = input('a for back, d for forward','s');
    nip = str2double(ip);
    if ip == 'a'
        image_pair_number = max(0,image_pair_number-1);
    elseif ip == 'd'
        image_pair_number = min(number_of_pair-1, image_pair_number+1);
    elseif ~isnan(nip) && floor(nip) == nip && image_pair_number >= 0 && image_pair_number < number_of_pair
        image_pair_number = nip;
    end

end



for image_pair_number = 123%0:number_of_pair-1%140%123%image_pair_number=0:number_of_pair-1
%PIV
temp = load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = temp.imgPiv;
temp = load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = temp.imgPiv;

[h, w] = size(IM_a);
 
%PIV Surf
surfOffset = 0;%10
imSurfa = FindSurfaceCapillary([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat'], findMask = true, surfOffset=surfOffset); 

imSurfb = FindSurfaceCapillary([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat'], findMask = true, surfOffset=surfOffset);

cl = [1,0.4,0.4];

figure(5)
hold off
imagesc(imSurfa.ImgScaledToPIVSmallCrop,[0,300])
hold on
plot(imSurfa.surfaceSurfImgScaled, '-', 'Color', cl)
plot(imSurfa.surface_raw, '-r')
daspect([1,1,1])

figure(6)
hold off
imagesc(imSurfb.ImgScaledToPIVSmallCrop,[0,300])
hold on
plot(imSurfb.surfaceSurfImgScaled, '-', 'Color', cl)
plot(imSurfb.surface_raw, '-r')
daspect([1,1,1])

figure(7)
imagesc(IM_a,[0,300])
% imagesc(IM_a.*imSurfa.mask,[0,300])
hold on
colormap gray
plot(imSurfa.surfacePIVImg,'Color', cl,'LineWidth',2)
daspect([1,1,1])
xlim([500,900])
ylim([200,600])

figure(8)
imagesc(IM_b,[0,300])
% imagesc(IM_b.*imSurfb.mask,[0,300])
hold on
colormap gray
plot(imSurfb.surfacePIVImg, 'Color', cl,'LineWidth',2)
daspect([1,1,1])
xlim([500,900])
ylim([200,600])

% imwrite(uint8(IM_a.*imSurfa.mask), sprintf('imgs/%s_%d_imgA.tif',exp_name,image_pair_number))
% imwrite(uint8(IM_b.*imSurfb.mask), sprintf('imgs/%s_%d_imgB.tif',exp_name,image_pair_number))

%% Compute Velocities
%IntrWndw=[256 128 64 32 16 8];
%GrdSpc=[128 64 32 16 8 4];
IntrWndw=[256 128 64 24 16 8];
GrdSpc=   [128 64 32 12 8 4];
% compVel = computeVelocities_marc_quick_nofilt(IM_a, IM_b, mask1, mask2, IntrWndw, GrdSpc);
compVel = ComputeVelocities_Quick_NoFilt_Deform_Water(IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc);
compVel.DX=1/17697.69; %m per pix
compVel.DT=10d-3; % sec per image pair -  DELTA_T= 10 milisec
%% Plot Deformed Image Pair
% [X,Y] = meshgrid((1:size(compVel.delx,2))*4,(1:size(compVel.delx,1))*4');
[X1,Y1] = meshgrid(1:w, 1:h);
% U1 = interp2(X,Y,compVel.delx,X1,Y1,'*spline'); 
% V1 = interp2(X,Y,compVel.dely,X1,Y1,'*spline');
% IM_b_D = interp2(1:size(IM_b,2),(1:size(IM_b,1))',IM_b,X1+U1,Y1-V1,'*linear');
U1 = compVel.delta_x1;
V1 = -compVel.delta_z1;
IM_b_D = interp2(1:size(IM_b,2),(1:size(IM_b,1))',IM_b,X1+U1,Y1+V1,'*linear');
figure(9)
imagesc(IM_b_D,[0,300])
colormap gray
daspect([1,1,1])
xlim([500,900])
ylim([200,600])
%% Plot PIV
[h,w] = size(compVel.delta_x);

figure(3)
hold off
% imagesc((1:w)*compVel.DX, (1:h)*compVel.DX, compVel.delta_x.*compVel.mask*compVel.DX/compVel.DT)
imagesc(compVel.xPIV*compVel.DX, compVel.zPIV*compVel.DX, compVel.delta_x.*compVel.mask*compVel.DX/compVel.DT)
ybot = 200;
% imagesc((1:w)*compVel.DX, (1:ybot*4)*compVel.DX, compVel.delx(1:ybot,:).*compVel.mask(1:ybot,:)*compVel.DX/compVel.DT)

hold on
colormap gray
set(gca, 'TickLabelInterpreter','latex', 'FontSize', 20);
c = colorbar;
c.TickLabelInterpreter = 'Latex';
c.Label.String = 'u (m/s)';
c.Label.Interpreter = 'latex';
xlabel('x (m)', 'Interpreter','latex')
ylabel('z (m)', 'Interpreter', 'latex')
% titleStr = sprintf('ExpLC%s-%s Pair Number %d', exp_name(end-3:end-3), exp_name(end-1:end), image_pair_number);
% title(titleStr,'Interpreter','latex')
daspect([1,1,1])

figure(4)
hold off
imagesc((1:w)*compVel.DX, (1:h)*compVel.DX, compVel.delta_z.*compVel.mask*compVel.DX/compVel.DT)
hold on
colormap gray
set(gca, 'TickLabelInterpreter','latex', 'FontSize', 20);
c = colorbar;
c.TickLabelInterpreter = 'Latex';
c.Label.String = 'w (m/s)';
c.Label.Interpreter = 'latex';
xlabel('x (m)', 'Interpreter','latex')
ylabel('z (m)', 'Interpreter', 'latex')
% titleStr = sprintf('ExpLC%s-%s Pair Number %d', exp_name(end-3:end-3), exp_name(end-1:end), image_pair_number);
% title(titleStr,'Interpreter','latex')
daspect([1,1,1])
%% Draw lines of constant s and n on the velocity field
surf = imSurfa.surfacePIVImg;
surfLen = length(surf);
s = 1:surfLen;
T = zeros(2,surfLen);
N = zeros(2,surfLen);
n = 0:20:500;
[ss,nn] = meshgrid(s,n);

T(:,1) = [1; surf(2)-surf(1)]/norm([1; surf(2)-surf(1)]);
N(:,1) = [-surf(2)+surf(1),1]/norm([-surf(2)+surf(1),1]);
for i = 2:surfLen-1
    T(:,i) = [1;(surf(i+1)-surf(i-1))/2]/norm([1;(surf(i+1)-surf(i-1))/2]);
    N(:,i) = [-(surf(i+1)-surf(i-1))/2;1]/norm([1;(surf(i+1)-surf(i-1))/2]);
end
T(:,end) = [1; surf(end)-surf(end-1)]/norm([1; surf(end)-surf(end-1)]);
N(:,end) = [-surf(end)+surf(end-1),1]/norm([-surf(end)+surf(end-1),1]);

[xsngrid,ysngrid,N] = sntoxy(ss,nn,surf);
figure(3)
for i = 1:length(n)
    plot(xsngrid(i,:)*compVel.DX,ysngrid(i,:)*compVel.DX,'-r')
    hold on
    % set(gca, 'YDir', 'reverse');
end
for i = 1:20:length(s)
    plot(xsngrid(:,i)*compVel.DX,ysngrid(:,i)*compVel.DX,'-r')
    hold on
end
%% Calculate xy position of each s,n point image A and image B and plot distorted images
surfa = imSurfa.surfacePIVImg;
surfb = imSurfb.surfacePIVImg;
s = 1:length(surfa);
n = 0:1:80;
[ss,nn] = meshgrid(s,n);

[xsngrida,ysngrida] = sntoxy(ss,nn,surfa);
[xsngridb,ysngridb] = sntoxy(ss,nn,surfb);
figure(11)
hold off
imagesc(IM_a,[0,300])
colormap gray
hold on
daspect([1,1,1])
for i = 1:length(n)
    plot(xsngrida(i,:),ysngrida(i,:),'-r')
    hold on
    % set(gca, 'YDir', 'reverse');
end
for i = 1:length(s)
    plot(xsngrida(:,i),ysngrida(:,i),'-r')
    hold on
end
IM_a_sn = interp2(IM_a,xsngrida,ysngrida);
fa = figure('Position',[1200,1200,1200,1200])
hold off
imagesc(s*compVel.DX, n*compVel.DX,IM_a_sn,[0,300])
hold on
xlim([0.01,0.03])
ylim([0,0.008])
colormap gray
daspect([1,1,1])

IM_b_sn = interp2(IM_b,xsngridb,ysngridb);
fb = figure('Position',[1200,1200,1200,1200])
hold off
imagesc(s*compVel.DX, n*compVel.DX, IM_b_sn,[0,300])
hold on
xlim([0.01,0.03])
ylim([0,0.008])
colormap gray
daspect([1,1,1])

% Converting u and v to sdot and ndot. Not because it is useful but so we
% can demonstrate how to go back the other way. We know x at a bunch of
% sn points, but we want to know s at a bunch of x,y points. 
% First, write the values of s and n on an xy grid
x = s;
y = 350:(350+80);
[X,Y] = meshgrid(x,y);
% At time A
sxygridaSI = scatteredInterpolant(xsngrida(:),ysngrida(:),ss(:));
sxygrida = sxygridaSI(X,Y);
nxygridaSI = scatteredInterpolant(xsngrida(:),ysngrida(:),nn(:));
nxygrida = nxygridaSI(X,Y);
% At time B
sxygridbSI = scatteredInterpolant(xsngridb(:),ysngridb(:),ss(:));
sxygridb = sxygridbSI(X,Y);
nxygridbSI = scatteredInterpolant(xsngridb(:),ysngridb(:),nn(:));
nxygridb = nxygridbSI(X,Y);

% Find partial derivatives of s with respect to x, y, and t, if s is
% defined as a function of these variables
[pspxa,pspya] = gradient(sxygrida);
[pnpxa,pnpya] = gradient(nxygrida);
pspt = sxygridb-sxygrida;
pnpt = nxygridb-nxygrida;

% Calculate sdot and ndot on an xy grid based on u and v on that grid
sdotxygrid =    pspxa.*compVel.delta_x1(y(1):y(end),x(1):x(end)) + ...
                pspya.*(-compVel.delta_z1(y(1):y(end),x(1):x(end))) + ...
                pspt;
ndotxygrid =    pnpxa.*compVel.delta_x1(y(1):y(end),x(1):x(end)) + ...
                pnpya.*(-compVel.delta_z1(y(1):y(end),x(1):x(end))) + ...
                pnpt;
% Convert sdot and ndot to be on an sn grid
sdotsngridSI = scatteredInterpolant(sxygrida(:),nxygrida(:),sdotxygrid(:));
sdotsngrid = sdotsngridSI(ss,nn);
ndotsngridSI = scatteredInterpolant(sxygrida(:),nxygrida(:),ndotxygrid(:));
ndotsngrid = ndotsngridSI(ss,nn);

% Now convert ndot and sdot back to u and v. We already know the values of
% x and y on the sn grid (xsngrida, ysngrida, etc.), so we can skip that
% step going back the other way. Now we need to find the partial
% derivatives of x and y with respect to s, n, and t, if x and y are
% defined as functions of these variables.
[pxpsa, pxpna] = gradient(xsngrida);
[pypsa, pypna] = gradient(ysngrida);
pxpt = xsngridb - xsngrida;
pypt = ysngridb - ysngrida;

% Calculated u and v on an sn grid based on sdot and ndot on that grid
usngrid = pxpsa.*sdotsngrid + pxpna.*ndotsngrid + pxpt;
vsngrid = pypsa.*sdotsngrid + pypna.*ndotsngrid + pypt;

% Convert u and v to be on an xy grid
uxygridSI = scatteredInterpolant(xsngrida(:),ysngrida(:),usngrid(:));
uxygrid = uxygridSI(X,Y);
vxygridSI = scatteredInterpolant(xsngrida(:),ysngrida(:),vsngrid(:));
vxygrid = vxygridSI(X,Y);

%% Try again on the sn transformation but starting from a synthetic sdot and ndot for surface gravity waves derived in Mathematica nb.
lambda = 100;
k = 2*pi/lambda;
ep = 1;
g = 1;
s = 1:(4*lambda);
n = fliplr(-lambda:0);
[ss,nn] = meshgrid(s,n);
dt = 1;
ndotsngrid = (-(-1 + exp(k*nn)).*sqrt(g*k).*sin(-k*ss)*ep - 0.5*ep^2*k*sqrt(g*k)*(-2 + exp(k*nn).*(2+k*nn)).*sin(-2*k*ss));
sdotsngrid = (sqrt(g*k) * (exp(k*nn) + k*nn).*cos(-k*ss)*ep - 0.5*ep^2*k*sqrt(g*k)*(-1 + k*nn.*(2*exp(k*nn)+k*nn)+(1-2*exp(k*nn)+k*nn.*(-4+k*nn)).*cos(-2*k*ss)));
figure(1)
hold off
imagesc(ndotsngrid)
% set(gca,'YDir','normal')
daspect([1,1,1])
figure(2)
hold off
imagesc(sdotsngrid)
daspect([1,1,1])


surfa = ep*cos(-k*s) + 0.5*k*ep^2*cos(-2*k*s);
surfb = ep*cos(sqrt(g*k)*dt-k*s) + 0.5*k*ep^2*cos(2*sqrt(g*k)*dt-2*k*s);
figure(3)
hold off
plot(surfa)
hold on
plot(surfb)

[xsngrida,ysngrida] = sntoxy(ss,nn,surfa);
[xsngridb,ysngridb] = sntoxy(ss,nn,surfb);
figure(11)
hold off

for i = 1:length(n)
    plot(xsngrida(i,:),ysngrida(i,:),'-r')
    hold on
    % set(gca, 'YDir', 'reverse');
end
for i = 1:length(s)
    plot(xsngrida(:,i),ysngrida(:,i),'-r')
    hold on
end
% imagesc(IM_a,[0,300])
colormap gray
daspect([1,1,1])

[~,gradn] = gradient(nn);
[grads,~] = gradient(ss);
[pxpsa, pxpna] = gradient(xsngrida);
pxpsa = pxpsa./grads; pxpna = pxpna./gradn;
[pypsa, pypna] = gradient(ysngrida);
pypsa = pypsa./grads; pypna = pypna./gradn;
pxpt = (xsngridb - xsngrida)/dt;
pypt = (ysngridb - ysngrida)/dt;

% Calculated u and v on an sn grid based on sdot and ndot on that grid
usngrid = pxpsa.*sdotsngrid + pxpna.*ndotsngrid + pxpt;
vsngrid = pypsa.*sdotsngrid + pypna.*ndotsngrid + pypt;

% Convert u and v to be on an xy grid
x = s;
y = fliplr((-lambda+ep+1):(ep));
[X,Y] = meshgrid(x,y);
uxygridSI = scatteredInterpolant(xsngrida(:),ysngrida(:),usngrid(:));
uxygrid = uxygridSI(X,Y);
vxygridSI = scatteredInterpolant(xsngrida(:),ysngrida(:),vsngrid(:));
vxygrid = vxygridSI(X,Y);

figure(4)
hold off
imagesc(uxygrid)
daspect([1,1,1])
figure(5)
hold off
imagesc(vxygrid)
daspect([1,1,1])
%% Save Results
% resultsPath = [load_path, '/Results_Surflab/'];
% tifPath = [resultsPath, 'tif/'];
% if ~exist(tifPath, 'dir')
%     mkdir(tifPath);
% end
% 
% compVel = rmfield(compVel, {'Mask', 'delx', 'dely', 'INTdelx', 'INTdelz'})
% imSurfaSave = rmfield(imSurfa, {'ImgScaledCroppedToPIV', 'ImgScaledToPIVSmallCrop', 'ImgScaledToPIV','mask'})
% imSurfbSave = rmfield(imSurfb, {'ImgScaledCroppedToPIV', 'ImgScaledToPIVSmallCrop', 'ImgScaledToPIV','mask'})
% 
% results_struct = struct();
% results_struct.compVel = compVel;
% results_struct.imSurfa = imSurfaSave;
% results_struct.imSurfb = imSurfbSave;
% 
% outputName = sprintf('%s/%s_%d_PIV.mat', resultsPath, exp_name, image_pair_number);
% save(outputName,'-v7.3','-fromstruct', results_struct)
% 
% imwrite(uint8(IM_a.*imSurfa.mask), sprintf('%s%s_%d_imgA.tif', tifPath, exp_name, image_pair_number))
% imwrite(uint8(IM_b.*imSurfb.mask), sprintf('%s%s_%d_imgB.tif', tifPath, exp_name, image_pair_number))
% 
% if ~exist(tifPath, 'dir')
%     mkdir(tifPath);
% end
% 
% y_crop = min([imSurfa.surface, imSurfb.surface]);
% 
% imwrite(uint8(IM_a),[tifPath, exp_name,'-',int2str(image_pair_number), '_a.tif'])
% imwrite(uint8(IM_b),[tifPath, exp_name,'-',int2str(image_pair_number), '_b.tif'])
% imwrite(uint8(IM_a.*maska),[tifPath, exp_name,'-',int2str(image_pair_number), '_masked_a.tif'])
% imwrite(uint8(IM_b.*maskb),[tifPath, exp_name,'-',int2str(image_pair_number), '_masked_b.tif'])
% imwrite(uint8(IM_a(y_crop:end,:).*maska(y_crop:end,:)),[tifPath, exp_name,'-',int2str(image_pair_number), '_masked_cropped_a.tif'])
% imwrite(uint8(IM_b(y_crop:end,:).*maskb(y_crop:end,:)),[tifPath, exp_name,'-',int2str(image_pair_number), '_masked_cropped_b.tif'])
% imwrite(uint8(IM_a),[tifPath, exp_name,'-',int2str(image_pair_number), '_shifted_a.tif'])
% imwrite(uint8(IM_b_D),[tifPath, exp_name,'-',int2str(image_pair_number), '_shifted_b.tif'])
% 
% %% Plot MLPIV velocities
% outfname = [tifPath, exp_name,'-',int2str(image_pair_number),'_shifted_a_out.flo']; %For Pairnum 0399
% [u, v] = read_flo_file(outfname);
% figure(10)
% hold off
% imagesc(compVel.delta_z1.*maska*compVel.DX/compVel.DT+u*compVel.DX/compVel.DT);
% hold on
% colormap gray
% set(gca, 'TickLabelInterpreter','latex', 'FontSize', 20);
% c = colorbar;
% c.TickLabelInterpreter = 'Latex';
% c.Label.String = 'w (m/s)';
% c.Label.Interpreter = 'latex';
% xlabel('x (m)', 'Interpreter','latex')
% ylabel('z (m)', 'Interpreter', 'latex')
% titleStr = sprintf('ExpLC%s-%s Pair Number %d', exp_name(end-3:end-3), exp_name(end-1:end), image_pair_number);
% title(titleStr,'Interpreter','latex')
% daspect([1,1,1])
end
end
%% Save results

% outfile = [load_path '/PIVMat_2023/' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
% save(outfile, 'compVel', 'imSurf1', 'imSurf2');
% disp(['pair ' num2str(image_pair_number) ' velocity done.']);

%%


% imSurf = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
% Surface_PIV=imSurf.surface;
% 

Surface_PIV = imSurfa.surfacePIVImg;
pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;
% 
% imwrite(uint8(IM_a.*imSurfa.mask), sprintf('imgs/%s_%d_imgA.tif',exp_name,image_pair_number))
% imwrite(uint8(IM_b.*imSurfb.mask), sprintf('imgs/%s_%d_imgB.tif',exp_name,image_pair_number))

transfo = generateTransfo_LC_noLFV_2023( compVel, Surface_PIV, pivRes); % 0 is there to compare SU(0,:) with surface
% 
SU = transfo.SU; 
SU = SU(2:end,:); % all but surface; % first line is zeta=0 THE SURFACE EXACTELY
% SU = SU -1716+287;
ORBX = transfo.ORBX;
ORBX = ORBX(2:end,:);
ORBZ = transfo.ORBZ;
ORBZ = ORBZ(2:end,:);
%
pivRes.GS = compVel.GS;
pivRes.zPIV = compVel.zPIV;
pivRes.pf_surf =SU(1,:);

u = compVel.delx.*compVel.mask;
w = compVel.dely.*compVel.mask;
%
intrp_u = transformVelField_decay_forFab( u, pivRes, SU );
intrp_w = transformVelField_decay_forFab( w, pivRes, SU );
% 
% intU_minusORBX = intrp_u - ORBX;
% intW_minusORBZ = intrp_w - ORBZ;
% 
% u_turb = reverseTransformVelField_decay_forFab( intU_minusORBX, pivRes, SU );
% w_turb = reverseTransformVelField_decay_forFab( intW_minusORBZ, pivRes, SU );
% uuTest = reverseTransformVelField_decay_forFab( intrp_u, pivRes, SU);


%% Function to convert from sn coordinates to xy coordinates
% assuming that s = x on the surface, and N is the normal vectors for the
% s coordinates given, and surf is the surface elevation for the given s
% coordinates.

function [x,y,N] = sntoxy(s,n,surf)
    surfLen = length(surf);
    T = zeros(2,surfLen);
    N = zeros(2,surfLen);
    
    T(:,1) = [1; surf(2)-surf(1)]/norm([1; surf(2)-surf(1)]);
    N(:,1) = [-surf(2)+surf(1),1]/norm([-surf(2)+surf(1),1]);
    for i = 2:surfLen-1
        T(:,i) = [1;(surf(i+1)-surf(i-1))/2]/norm([1;(surf(i+1)-surf(i-1))/2]);
        N(:,i) = [-(surf(i+1)-surf(i-1))/2;1]/norm([1;(surf(i+1)-surf(i-1))/2]);
    end
    T(:,end) = [1; surf(end)-surf(end-1)]/norm([1; surf(end)-surf(end-1)]);
    N(:,end) = [-surf(end)+surf(end-1),1]/norm([-surf(end)+surf(end-1),1]);
    
    x = zeros(size(s));
    y = zeros(size(n));
    
    for i = 1:size(n,1)
        x(i,:) = s(i,:) + n(i,:).*N(1,:);
        y(i,:) = surf + n(i,:).*N(2,:);
    end
end

function [modstruct] = rmfields(ogstruct, fields_to_delete)
    existing_fields = intersect(fields_to_delete, ...
        fieldnames(ogstruct));
    if ~isempty(existing_fields)
        ogstruct = rmfield(ogstruct, existing_fields);
        fprintf('  Removed %d field(s)', ...
            numel(existing_fields));
    else
        fprintf('  No matching fields to remove from %s in %s\n', ...
            sname, fname);
    end
    modstruct = ogstruct;
end