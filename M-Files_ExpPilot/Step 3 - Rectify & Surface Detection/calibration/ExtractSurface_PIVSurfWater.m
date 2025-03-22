function [BadFramePIVSurfW,XPIVSurfW_Surface,PIVSurfW_Surface,XPIVW_PIVSurfW_Surface,PIVW_PIVSurfW_Surface,PIVW_Surface] = ExtractSurface_PIVSurfWater(PIVSurfW_CamAngle,PIV_W)

%% PIVSurf Water surface detection

%% Step 1: Find surface
% Extract surface
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

%% Step 2 : Resize PIVSurf image to the Size of PIVWater images
%%% Find landmarks of PIV Water in PIVSurf Water coordinates 
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

T2 = maketform('projective',U1,X1);

[Xsurf2,Ysurf2] = tformfwd(T2,Usurf,Vsurf);

% Further rototranslation to match flat surface
RotAngle = -39/6107;
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf2 = M*[Xsurf2;Ysurf2];
Ysurf2 = Ysurf2+17;

Xsurf = round(Ysurf2(1,1)):round(Ysurf2(1,end));
Ysurf = interp1(Xsurf2,Ysurf2(2,:),Xsurf,'linear','extrap');

%% Step 3: Save Surfaces in all coordinate systems
%%% PIVSurf Water Surface in PIVSurf Water coordinates 
XPIVSurfW_Surface = Usurf;
PIVSurfW_Surface = Vsurf;

%%% PIVSurf Water Surface in PIVWater coordinates
XPIVW_PIVSurfW_Surface = Xsurf;
PIVW_PIVSurfW_Surface = Ysurf;

[~,Ix] = min(abs(Xsurf-4));
PIVW_Surface = PIVW_PIVSurfW_Surface(Ix:size(PIV_W,2)+Ix-1);
