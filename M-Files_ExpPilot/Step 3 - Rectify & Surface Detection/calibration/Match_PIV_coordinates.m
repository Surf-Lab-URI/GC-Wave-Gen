function [CompVelWater_PIVAir,CompVelWater_PIVSurfW,Uinv,Vinv] = Match_PIV_coordinates(XPIVW_PIVSurfW1_Surface,PIVW_PIVSurfW1_Surface,XPIVSurfW1_Surface,PIV_LFV_Surface,CompVelWater,PIV1_W,PIV1_A)

%% Transform PIV Water measurements in PIVSurf Water coordinates

%% Step 1: Inverse transform of the SURFACE - From PIV Water to PIVSurf Water
% This step is necessary to verify the matching between original surface and the transformed one 

%%% Transformation from PIV Water to PIVSurf Water coordinate system of the
%%% surface
% Points to use for projection of Landmarks on the PIVSurf Water grid
Rp1 = fliplr([ 327 3058 ; 731 3056 ; 1135 3052 ; 1541 3048 ; 1944 3045 ]); % GRID POINTS on the first straight line RIGHT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_R = -3.25; DY_R = 404.25; % main DIFFERENCE in the horizontal (DX_R) and in the vertical (DY_R) 
Lp1 = [1013 311 ; 1010 714 ; 1007 1118 ; 1004 1523 ; 1000 1927 ]; % GRID POINTS on the first straight line LEFT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_L = -3.25; DY_L = 404; % main DIFFERENCE in the horizontal (DX_L) and in the vertical (DY_L) 
% Landmarks in the PIVSurf Water grid 
Xp1(2:3,1:2) = [Rp1(end,1)+DX_R Rp1(end,2)+DY_R ; Rp1(end,1)+4*DX_R Rp1(end,2)+4*DY_R];
Xp1([1,4],1:2) = [Lp1(end,1)+DX_L Lp1(end,2)+DY_L ; Lp1(end,1)+4*DX_L Lp1(end,2)+4*DY_L];

U1 = Xp1; % COORDINATES IN PIVSURF WATER
X1 = [ 404 1195 ; 3447 1198 ; 3447 3018 ; 405 3019 ]; % COORDINATES IN PIV WATER

T1 = maketform('projective',U1,X1); %Transformation from PIVSurf Water to PIV Water

% Eliminate the rototranslation added after the transformation
DY = 17; % Same translation used in Extract_PIVSurf Water
RotAngle = -39/6107; % Same rotation angle used in Extract_PIVSurf Water
Minv = [cos(-RotAngle) sin(-RotAngle); -sin(-RotAngle) cos(-RotAngle)];
Yinv2 = PIVW_PIVSurfW1_Surface-DY;
Yinv2 = Minv*[XPIVW_PIVSurfW1_Surface;Yinv2];
Xinv = XPIVW_PIVSurfW1_Surface;
Yinv = Yinv2(2,:);

% Inverse mapping of the surface
[Uinv2,Vinv2] = tforminv(T1,Xinv,Yinv);
Uinv = XPIVSurfW1_Surface;
Vinv = interp1(Uinv2,Vinv2,Uinv,'linear','extrap');

%%% Original surface: XPIVSurfW1_Surface,PIVSurfW1_Surface
%%% Transformed surface: Uinv, Vinv
% figure;imagesc(PIVSurfW1_CamAngle);colormap gray;caxis([0 200]);hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')
% figure;plot(Uinv,PIVSurfW1_Surface-Vinv,'r')

%% Step 2: Inverse transform of the PIV Water VELOCITY - From PIV Water to PIVSurf Water

%%% Inverse transformation (from PIV Water to PIVSurf Water)
% Points to use for projection of Landmarks on the PIVSurf Water grid
Rp1 = fliplr([ 327 3058 ; 731 3056 ; 1135 3052 ; 1541 3048 ; 1944 3045 ]); % GRID POINTS on the first straight line RIGHT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_R = -3.25; DY_R = 404.25; % main DIFFERENCE in the horizontal (DX_R) and in the vertical (DY_R) 
Lp1 = [1013 311 ; 1010 714 ; 1007 1118 ; 1004 1523 ; 1000 1927 ]; % GRID POINTS on the first straight line LEFT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_L = -3.25; DY_L = 404; % main DIFFERENCE in the horizontal (DX_L) and in the vertical (DY_L) 
% Landmarks in the PIVSurf Water grid 
Xp1(2:3,1:2) = [Rp1(end,1)+DX_R Rp1(end,2)+DY_R ; Rp1(end,1)+4*DX_R Rp1(end,2)+4*DY_R];
Xp1([1,4],1:2) = [Lp1(end,1)+DX_L Lp1(end,2)+DY_L ; Lp1(end,1)+4*DX_L Lp1(end,2)+4*DY_L];

U1 = Xp1; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [ 404 1195 ; 3447 1198 ; 3447 3018 ; 405 3019 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.

T2inv = maketform('projective',X1,U1);

[CompVelWater_PIVSurfW] = transform_PIVWater_to_PIVSurf_Water(CompVelWater,PIV1_W,T2inv);

%%% Transformed PIV Water velocity: CompVelW_inv
% figure;imagesc(XposInv(1):XposInv(end),YposInv(1):YposInv(end),CompVelW_inv);
% hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')

%% Step 3: Transform surface and PIV Water velocity - from PIVSurf Water to PIV Air coordinates
U2 = [999 1927 ; 1009 816 ; 3056 832 ; 3045 1944 ]; % The coordinate
% of four points in the PIVSurf Air - LFV surface image (the grid calibration image).
X2 = [481 1877 ; 490 44 ; 3825 63 ; 3815 1894 ]; % The coordinates of  the four
% points in the PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T3 = maketform('projective',U2,X2);

[CompVelWater_PIVAir] = transform_PIVSurf_Water_to_PIVAir(T3,Uinv,Vinv,CompVelWater_PIVSurfW,PIV_LFV_Surface,PIV1_A);
