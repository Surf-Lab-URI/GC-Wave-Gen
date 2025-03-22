function [CompVelWater_PIVAir,CompVelWater_PIVSurfW,Uinv,Vinv] = Match_PIV_coordinates(RotAngle_W,DY_W,T2_W,XPIVW_PIVSurfW1_Surface,PIVW_PIVSurfW1_Surface,XPIVSurfW1_Surface,XPIV_LFV_Surface,PIV_LFV_Surface,CompVelWater,PIV1_W,PIV1_A)

%% Transform PIV Water measurements in PIVSurf Water coordinates

%%% Inverse transformation (from PIV Water to PIVSurf Water)
% Points to use for projection of Landmarks on the PIVSurf Water grid
%%% Projection of the grid points under the water
% Left points LP1 and LP2
LP_Inch8th = [1203 1963 ; 1203 1915 ; 1204 1866 ; 1204 1819 ; 1204 1772 ; 1205 1724 ; 1205 1675 ; 1206 1627 ; 1206 1579 ; 1207 1531 ; 1207 1483 ; 1207 1436 ; 1208 1387 ; 1208 1338 ; 1208 1290 ];
DeltaY_LP = -mean(diff(LP_Inch8th(:,2)));
CC1 = polyfit(LP_Inch8th(:,2),LP_Inch8th(:,1),1); % straight line fitting landmarks
LP1(1,2) = LP_Inch8th(1,2)+4*DeltaY_LP;
LP1(1,1) = polyval(CC1,LP1(1,2));
LP2(1,2) = LP1(1,2)+(3*8)*DeltaY_LP;
LP2(1,1) = polyval(CC1,LP2(1,2));
% Right points RP1 and RP2
RP_Inch8th = [ 3163 1986 ; 3163 1938 ; 3164 1890 ; 3164 1841 ; 3164 1793 ; 3165 1745 ; 3165 1697 ; 3165 1649 ; 3166 1601 ; 3166 1553 ; 3166 1503 ; 3167 1455 ;  3167 1407 ; 3168 1359 ;  3168 1311];
DeltaY_RP = -mean(diff(RP_Inch8th(:,2)));
CC2 = polyfit(RP_Inch8th(:,2),RP_Inch8th(:,1),1); % straight line fitting landmarks
RP1(1,2) = RP_Inch8th(1,2)+4*DeltaY_RP;
RP1(1,1) = polyval(CC2,RP1(1,2));
RP2(1,2) = RP1(1,2)+(3*8)*DeltaY_RP;
RP2(1,1) = polyval(CC2,RP2(1,2));


U2 = [ LP1 ; LP2 ; RP2 ; RP1]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X2 = [ 630 875 ; 610 2715 ; 3706 2743; 3718 903 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.

%% Step 1: Inverse transform of the SURFACE - From PIV Water to PIVSurf Water
% This step is necessary to verify the matching between original surface and the transformed one 

%%% Transformation from PIV Water to PIVSurf Water coordinate system of the
%%% surface
% Eliminate the rototranslation added after the transformation
Minv = [cos(-RotAngle_W) sin(-RotAngle_W); -sin(-RotAngle_W) cos(-RotAngle_W)];
Yinv2 = PIVW_PIVSurfW1_Surface+DY_W;
Yinv2 = Minv*[XPIVW_PIVSurfW1_Surface;Yinv2];
Xinv = Yinv2(1,:);
Yinv = Yinv2(2,:);
I = find(isnan(Xinv));
% Delete NaN coming out after rotation
Xinv(I) = [];
Yinv(I) = [];

% T2 = maketform('projective',U2,X2);
% T2 = fitgeotform2d(U2,X2,"projective");

% Inverse mapping of the surface
[Uinv2,Vinv2] = tforminv(T2_W,Xinv,Yinv);
Uinv = XPIVSurfW1_Surface;
Vinv = interp1(Uinv2,Vinv2,Uinv,'linear','extrap');

%%% Original surface: XPIVSurfW1_Surface,PIVSurfW1_Surface
%%% Transformed surface: Uinv, Vinv
% figure;imagesc(PIVSurfW1_CamAngle);colormap gray;caxis([0 200]);hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')
% figure;plot(Uinv,PIVSurfW1_Surface-Vinv,'r')

%% Step 2: Inverse transform of the PIV Water VELOCITY - From PIV Water to PIVSurf Water

T2inv = maketform('projective',X2,U2);

[CompVelWater_PIVSurfW] = transform_PIVWater_to_PIVSurf_Water(CompVelWater,PIV1_W,T2inv);

%%% Transformed PIV Water velocity: CompVelW_inv
% figure;imagesc(XposInv(1):XposInv(end),YposInv(1):YposInv(end),CompVelW_inv);
% hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')

%% Step 3: Transform surface and PIV Water velocity - from PIVSurf Water to PIV Air coordinates
U2 = [ 1203 1963 ; 1213 1001 ; 3171 1016; 3163 1986]; % The coordinate
% of four points in the PIVSurf Water surface image (the grid calibration image).
X2 = [ 595 2033 ; 615 454 ; 3780 480; 3768 2062 ]; % The coordinates of  the four
% points in the PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T3 = maketform('projective',U2,X2);

[CompVelWater_PIVAir] = transform_PIVSurf_Water_to_PIVAir(T3,Uinv,Vinv,CompVelWater_PIVSurfW,XPIV_LFV_Surface,PIV_LFV_Surface,PIV1_A);
