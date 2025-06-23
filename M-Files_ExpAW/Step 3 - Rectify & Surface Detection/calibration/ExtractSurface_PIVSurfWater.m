  function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface,XPIVW_PIVSurfW_Surface,PIVW_PIVSurfW_Surface,PIVW_Surface,T2,RotAngle,DY] = ExtractSurface_PIVSurfWater(PIVSurfW_CamAngle,PIV_W,Water_Surface)

%% PIVSurf Water surface detection

%% Step 1: Find surface
% Extract surface
YY = 501:size(PIVSurfW_CamAngle,2)-40;
XX = 1:2800;

%%% Pre-extrapolation
%PIVSurfW_CamAngle2 = adapthisteq(PIVSurfW_CamAngle/max(PIVSurfW_CamAngle(:)),'NumTiles',[8 8],'NBins',20); %reduce the sensitivity to eliminate noise

PIVSurfW_CamAngle2 = PIVSurfW_CamAngle;
DX = 200;
S_T = PIVSurfW_CamAngle2(XX(1+DX:end),YY);
S_B = PIVSurfW_CamAngle2(XX(1:end-DX),YY);
S2 = S_T-S_B;
% S = imfilter(S2,fspecial('gaussian',64,64),'replicate');

S = S2;
surfSigmas = [50,40,30,20,10,8,6];
surfSteps = [50,40,30,20,10,8];
surfMask = 1;
[imSurf] = CrapperOptimized_FindSurface(S,surfSigmas, surfSteps, surfMask); %FindSurface_Water(S, 5, 5);
PIVSurf_Surface_Raw = imSurf.surface;
% PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = PIVSurf_Surface_Raw;
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 3d2); %2d4 for looking at larger gravity waves
if length(SP.coefs)>2
    PIVSurf_Surface_W = SP.coefs(2:end-1);
else
    CC = polyfit([YY(1),YY(end)],[SP.coefs(1) SP.coefs(2)],1);
    PIVSurf_Surface_W = polyval(CC,YY(1):YY(end));
end
Usurf = YY;
Vsurf = PIVSurf_Surface_W+DX;

%%% Check if bad frame
BadFramePIVSurfW = 0;
if imSurf.badFrameBool == 1
    BadFramePIVSurfW = 1;
end

%% Step 2 : Resize PIVSurf image to the Size of PIVWater images
%%% Find landmarks of PIV Water in PIVSurf Water coordinates 
% Points used for PIV Water coordinates references
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

T2 = maketform('projective',U2,X2);

[Xsurf2,Ysurf2] = tformfwd(T2,Usurf,Vsurf);

% Further rototranslation to match flat surface
RotAngle = (674.955-684.981)/(5250.53+477.054);
DY = 4;%25;
if strcmp(Water_Surface,'2')
    DY = DY-8;
end
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf2 = M*[Xsurf2;Ysurf2];
Xsurf2 = Ysurf2(1,:);
Ysurf2 = Ysurf2(2,:);
Ysurf2 = Ysurf2-DY;

% Xsurf = round(Xsurf2(1):Xsurf2(end));
Xsurf = -477:5248; % hard-coded to make all the surfaces exactly the same length
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf,'spline','extrap');

%% Step 3: Save Surfaces in all coordinate systems
%%% PIVSurf Water Surface in PIVSurf Water coordinates 
XPIVSurfW_Surface = Usurf;
PIVSurfW_Surface = Vsurf;

%%% PIVSurf Water Surface in PIVWater coordinates
XPIVW_PIVSurfW_Surface = Xsurf;
PIVW_PIVSurfW_Surface = Ysurf;

[~,Ix] = min(abs(Xsurf-1));
PIVW_Surface = PIVW_PIVSurfW_Surface(Ix:size(PIV_W,2)+Ix-1);
