clear
LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1%:length(DIRS)

exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [LONG exp_name];
files=dir([load_path '/PIVRaw/PIV/*.mat']);
number_of_pair=length(files)/2;

%%%%%%%% Stuff added by Andy to quickly flip through frames. should be
%%%%%%%% commented out for normal use

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

%%%%%% End of stuff added by Andy to quickly flip through frames.

for image_pair_number = 105%123%image_pair_number=0:number_of_pair-1
%PIV
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imgPiv;
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = imgPiv;
 
%PIV Surf
load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurfa = imgPivsurf;
load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurfb = imgPivsurf;

%Surface detection and Creating Masks
U1 = [147 49;2024 57; 1995 1004; 161 999];
X1 = [147 49; 2024 49; 2024 1004; 147 1004];
T1 = fitgeotrans(U1,X1,'projective');
d1=imwarp(imgPivsurfa,T1,'cubic');
d2=imwarp(imgPivsurfb,T1,'cubic');

d1=imresize(d1,176.9769/105.5880);%Resizing to match PIV
d2=imresize(d2,176.9769/105.5880);%Resizing to match PIV
s1=d1(30:3525,755:755+2047); %croping
s2=d2(30:3525,755:755+2047);

%%%%%OG Fabrice surface detection
% imSurf1 = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
% imSurf2 = findSurface_simple_ext_force_2023((medfilt2(s2)), 1);

%%%%%Andy Crapper-optimized surface detection
surfSigmas = [50 40 30 20 15];
surfSteps = [50 40 30 5];
SurfMask = 1;
slopeDiffThreshold = 5;
imSurfa = CrapperOptimized_FindSurface(s1, surfSigmas, surfSteps, 1, slopeDiffThreshold);
imSurfa.surface = FiltSurf(imSurfa.surface_raw,200);
imSurfb = CrapperOptimized_FindSurface(s2, surfSigmas, surfSteps, 1, slopeDiffThreshold);
imSurfb.surface = FiltSurf(imSurfb.surface_raw,200);

cl = [1,0.4,0.4];

figure(5)
hold off
imagesc(s1,[0,300])
hold on
plot(imSurfa.surface, '-', 'Color', cl)
plot(imSurfa.surface_raw, '-r')
daspect([1,1,1])

figure(6)
hold off
imagesc(s2,[0,300])
hold on
plot(imSurfb.surface, '-', 'Color', cl)
plot(imSurfb.surface_raw, '-r')
daspect([1,1,1])

imSurfa.surface=imSurfa.surface-1716+287;
imSurfb.surface=imSurfb.surface-1716+287;
[h, w] = size(IM_a); %image height and width
maska=ones(size(IM_a));
maskb=ones(size(IM_a));
for i=1:w
    maska(1:round(imSurfa.surface(i)),i)=NaN;
    maskb(1:round(imSurfb.surface(i)),i)=NaN;
end
warning off
imSurfa.mask=maska;
imSurfb.mask=maskb;

figure(7)
imagesc(IM_a,[0,300])
hold on
colormap gray
plot(imSurfa.surface,'Color', cl,'LineWidth',2)
daspect([1,1,1])
xlim([500,900])
ylim([200,600])

figure(8)
imagesc(IM_b,[0,300])
hold on
colormap gray
plot(imSurfb.surface, 'Color', cl,'LineWidth',2)
daspect([1,1,1])
xlim([500,900])
ylim([200,600])

%% Compute Velocities
%IntrWndw=[256 128 64 32 16 8];
%GrdSpc=[128 64 32 16 8 4];
IntrWndw=[256 128 64 24 16 8];
GrdSpc=   [128 64 32 12 8 4];
% compVel = computeVelocities_marc_quick_nofilt(IM_a, IM_b, mask1, mask2, IntrWndw, GrdSpc);
compVel = ComputeVelocities_Quick_NoFilt_Deform_Water(IM_a, IM_b, maska, maskb, IntrWndw, GrdSpc);
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
figure(3)
hold off
imagesc((1:w)*compVel.DX, (1:h)*compVel.DX, compVel.delx.*compVel.mask*compVel.DX/compVel.DT)
hold on
colormap gray
set(gca, 'TickLabelInterpreter','latex', 'FontSize', 20);
c = colorbar;
c.TickLabelInterpreter = 'Latex';
c.Label.String = 'u (m/s)';
c.Label.Interpreter = 'latex';
xlabel('x (m)', 'Interpreter','latex')
ylabel('z (m)', 'Interpreter', 'latex')
titleStr = sprintf('ExpLC%s-%s Pair Number %d', exp_name(end-3:end-3), exp_name(end-1:end), image_pair_number);
title(titleStr,'Interpreter','latex')
daspect([1,1,1])

figure(4)
hold off
imagesc((1:w)*compVel.DX, (1:h)*compVel.DX, compVel.dely.*compVel.mask*compVel.DX/compVel.DT)
hold on
colormap gray
set(gca, 'TickLabelInterpreter','latex', 'FontSize', 20);
c = colorbar;
c.TickLabelInterpreter = 'Latex';
c.Label.String = 'w (m/s)';
c.Label.Interpreter = 'latex';
xlabel('x (m)', 'Interpreter','latex')
ylabel('z (m)', 'Interpreter', 'latex')
titleStr = sprintf('ExpLC%s-%s Pair Number %d', exp_name(end-3:end-3), exp_name(end-1:end), image_pair_number);
title(titleStr,'Interpreter','latex')
daspect([1,1,1])
%% Draw lines of constant s and n on the velocity field
surf = imSurfa.surface;
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

[xsngrid,ysngrid] = sntoxy(ss,nn,surf,N);
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
%% Generate tif files
% resultsPath = [load_path, '/Results_Surflab/'];
% tifPath = [resultsPath, '/MLPIVtif/'];
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
imSurf = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
Surface_PIV=imSurf.surface;

pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;

transfo = generateTransfo_LC_noLFV_2023( compVel, Surface_PIV, pivRes); % 0 is there to compare SU(0,:) with surface

SU = transfo.SU; 
SU = SU(2:end,:); % all but surface; % first line is zeta=0 THE SURFACE EXACTELY
SU = SU -1716+287;
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

intU_minusORBX = intrp_u - ORBX;
intW_minusORBZ = intrp_w - ORBZ;

u_turb = reverseTransformVelField_decay_forFab( intU_minusORBX, pivRes, SU );
w_turb = reverseTransformVelField_decay_forFab( intW_minusORBZ, pivRes, SU );
uuTest = reverseTransformVelField_decay_forFab( intrp_u, pivRes, SU);


%% Function to convert from sn coordinates to xy coordinates
% assuming that s = x on the surface, and N is the normal vectors for the
% s coordinates given, and surf is the surface elevation for the given s
% coordinates.

function [x,y] = sntoxy(s,n,surf,N)
x = zeros(size(s));
y = zeros(size(n));

for i = 1:size(n,1)
    x(i,:) = s(i,:) + n(i,:).*N(1,:);
    y(i,:) = surf + n(i,:).*N(2,:);
end
end


