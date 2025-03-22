function [PIVSurfW_CamAngle] = PIVSurfWater_CamAngle_Correction(IM)

%% PIVSurf Air - LFV CamAngle Correction

%% First: rotation to make surface flat
Angle = rad2deg(atan(6/4096));
PIVSurfW_rot = imrotate(IM,-Angle);

%% Second: affine transformation to rectify image
% Correct PIV Surf Camera Angle
U1 = [1052 2010 ; 1040 951 ; 3280.5 951 ; 3265 2010]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 =  [ (1052+1040)/2 2010 ; (1052+1040)/2 951 ; (3280.5+3265)/2 951 ; (3280.5+3265)/2 2010 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVSurfW_CamAngle =imwarp(PIVSurfW_rot,T1,'cubic');

%% Third: possible other rototranslation to match the images