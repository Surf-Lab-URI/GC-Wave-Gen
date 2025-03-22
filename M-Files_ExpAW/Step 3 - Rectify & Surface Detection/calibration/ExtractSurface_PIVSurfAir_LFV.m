function [BadFramePIVSurfLFV,XLFV_Surface,LFV_Surface,XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface] = ExtractSurface_PIVSurfAir_LFV(PIVSurfA_CamAngle,PIV1_A)

%% PIVSurf Air - LFV CamAngle Correction

%% Step 1: Find surface
% Extract surface
XX = 1001:3000;
YY = 251:8300;
% [imSurf] = Copy_of_FindSurface(PIVSurfA_CamAngle(1001:3000,:), 5, 5);
[imSurf] = Copy_of_FindSurface(PIVSurfA_CamAngle(XX,YY), 5, 5);
PIVSurf_Surface_Raw = imSurf.surface;
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1d4); %10000);
PIVSurf_Surface_A = SP.coefs(2:end-1);
% Usurf = 501:length(PIVSurf_Surface_A)-500;
% Vsurf = PIVSurf_Surface_A(501:end-500)+1000;
Usurf = YY;
Vsurf = PIVSurf_Surface_A+XX(1)-1;

%% Step 2 : Resize PIV Surf images to the Size of Fused PIV images
U2 = [3169 2384 ; 3169 1296 ; 5213 1299 ; 5210 2387]; % The coordinate
% of four points in the PIVSurf Air - LFV surface image (the grid calibration image).
X2 = [478 1874 ; 487 43 ; 3822 58 ; 3814 1892 ]; % The coordinates of  the four
% points in the PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T2 = maketform('projective',U2,X2);

%%% This is the actual resized image; we retrieve it only for checking,
%%% otherwise we transform only the surface (the resized image is huge!!)
% [Resized_PIVSurf,XPos,YPos] =  imtransform(PIVSurfA_CamAngle,T2,'XYScale',1);
%%% XPos and YPos retrieved from resized transformation!!
XPos(1) = -4.745351888322090e+03;
XPos(2) = 9.561648111677910e+03;
YPos(1) = -2.137815623096811e+03;
YPos(2) = 8.690184376903190e+03;
[Xsurf2,Ysurf2] = tformfwd(T2,Usurf,Vsurf);

% Further rotation to match flat surface
RotAngle = atan(30/4110);
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf2 = M*[1:length(Ysurf2);Ysurf2];
Ysurf2 = Ysurf2(2,:);
Ysurf2 = Ysurf2+20;

Xsurf = round(Xsurf2(1):Xsurf2(end));
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf);

% XPos and YPos are a two-element, real vector that  together  specifiy the
% spatial location of the output image B in the 2D output space XY. The two
% elements of XPos and YPos give the x-coordinates and y-coordinates of the
% first and last columns of B, respectively.

%%% Check if bad frame
BadFramePIVSurfLFV = 0;
if imSurf.badFrameBool == 1
    BadFramePIVSurfLFV = 1;
end

%%% LFV Surface in LFV coordinates ( == to PIVSurf Air)
XLFV_Surface = Usurf;
LFV_Surface = Vsurf;

%%% LFV Surface in PIV coordinates
XPIV_LFV_Surface = Xsurf;
PIV_LFV_Surface = Ysurf;

%%% PIV Surface from LFV
[~,Ix] = min(abs(Xsurf-1));
[~,Iy] = min(abs((YPos(1):YPos(2))-1));
PIV_Surface = Ysurf(Ix:Ix+4140-1);
%%% Resized PIVSurfA_LFV matching PIV image
% PIV_PIVSurfA = PIVSurfA_LFV_Corrected.img(Iy:Iy+3090-1,Ix:Ix+4140-1);

%% Remind : PIVSurf Air == LFV !!