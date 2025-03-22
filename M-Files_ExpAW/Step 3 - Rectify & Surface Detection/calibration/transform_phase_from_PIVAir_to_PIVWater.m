function [XPIVW_LFV_Surface,PIVW_LFV_Surface] = transform_phase_from_PIVAir_to_PIVWater(XPIV_LFV_Surface,PIV_LFV_Surface)
%% TRANSFORM LFV SURFACE TO PIVWATER COORDINATES 

%% From PIV Air to PIVSurf Water
XX2 = [ 1203 1963 ; 1213 1001 ; 3171 1016; 3163 1986]; % The coordinate
% of four points in the PIVSurf Water surface image (the grid calibration image).
U2 = [ 595 2033 ; 615 454 ; 3780 480; 3768 2062 ]; % The coordinates of  the four
% points in the PIV Air image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.
T3 = maketform('projective',U2,XX2);
Uinv = XPIV_LFV_Surface(~isnan(PIV_LFV_Surface));
Vinv = PIV_LFV_Surface(~isnan(PIV_LFV_Surface));
[XPIVSurfW_LFV_Surface,PIVSurfW_LFV_Surface] = tformfwd(T3,Uinv,Vinv);

%% From PIVSurf Water to PIVWater
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
% four corner of a quadrilateral in PIV Water
X2 = [ 630 875 ; 610 2715 ; 3706 2743; 3718 903 ]; % The coordinates of 
% four corner of quadrilateral in PIVSurf Water

T2 = maketform('projective',U2,X2);

[XX1,YY1] = tformfwd(T2,XPIVSurfW_LFV_Surface,PIVSurfW_LFV_Surface);
RotAngle = (674.955-684.981)/(5250.53+477.054);
DY = 4;%25;
%if strcmp(Water_Surface,'2')
%    DY = DY-8;
%end
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
YY2 = M*[XX1;YY1];
XX2 = YY2(1,:);
YY2 = YY2(2,:);
YY2 = YY2-DY;

% Xsurf = round(Xsurf2(1):Xsurf2(end));
%Xsurf = -477:5248; % hard-coded to make all the surfaces exactly the same length
XPIVW_LFV_Surface = round(XX2(1)):round(XX2(end));
PIVW_LFV_Surface = interp1(XX2,YY2,XPIVW_LFV_Surface,'spline','extrap');

%%% Further roto-translation to match the surface
DY2 = -11.5;
RotAngle2 = 11/(5248+477);
M2 = [cos(RotAngle2) sin(RotAngle2); -sin(RotAngle2) cos(RotAngle2)];
YY22 = M2*[XPIVW_LFV_Surface;PIVW_LFV_Surface];
XX22 = YY22(1,:);
YY22 = YY22(2,:);
YY22 = YY22-DY2;

%% Results
% XPIVW_LFV_Surface = round(XX22(1)):round(XX22(end));
XPIVW_LFV_Surface = -2346:7651; % hard-coded to make all the surfaces exactly the same length
PIVW_LFV_Surface = interp1(XX22,YY22,XPIVW_LFV_Surface,'spline','extrap');