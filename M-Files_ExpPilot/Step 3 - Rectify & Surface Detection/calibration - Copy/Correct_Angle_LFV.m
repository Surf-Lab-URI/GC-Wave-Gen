function [LFV_Corrected,T2_bis] = Correct_Angle_LFV(PIV2)
% Correct PIV Surf Camera Angle

%% Correct residual distortion through polynomial
%%% No further corrections for LFV

%% Correction to get straight plumblines and flat surface through projective
CamAngle = PIV2;
CamAngle(isnan(CamAngle)) = 0;
%%% Points retrieved from plumblines and flat surface
UU1 = [ 184.5  1457 ; 3788.5 1497 ; 3991.5 100 ; 1.75 100];
XX1 = [ (184.5+1.75)/2 1477 ; (3788.5+3991.5)/2 1477 ; (3788.5+3991.5)/2 100 ; (184.5+1.75)/2 100 ];

%%% Projective transformation to correct PIVSURF Camera Angle
%T1 = maketform('projective',U1,X1); % Creates spatial transformation struct
TT1 = fitgeotrans(UU1,XX1,'projective');
% for a two-dimensional projective transformation that map each row of U to
% the corresponding  row of X. The U and X arguments are each 4-by-2 matrix
% and define the corners of input and output quadrilaterals.  Note that  no
% three corners can be collinear.
PIV2_CamAngle_Corrected = imwarp(CamAngle,TT1,'linear'); %'cubic' instead of 'linear'?

%% Resize LFV Surf images to the Size of PIVSURF_CamAngle_Corrected.img

Up = [ 965 1026 ; 2947 876.5 ; 3459 1419 ; 968 1417 ];
Xp = [ 1196 1887 ; 7951 1351 ; 9759 3182 ; 1190 3187 ];

T2_bis = maketform('projective',Up,Xp);
[Resized_LFV_bis,XPos,YPos] =  imtransform(PIV2_CamAngle_Corrected,T2_bis,'XYScale',1);

Resized_LFV = Resized_LFV_bis;

LFV_Corrected.img  = Resized_LFV;
LFV_Corrected.XPos = XPos;
LFV_Corrected.YPos = YPos;

% XPos and YPos are a two-element, real vector that  together  specifiy the
% spatial location of the output image B in the 2D output space XY. The two
% elements of XPos and YPos give the x-coordinates and y-coordinates of the
% first and last columns of B, respectively.

end