function [PIVWater_CamAngle] = PIVWater_CamAngle_Correction(IM)

%% PIV Air CamAngle Correction

%% First: rotation to make surface flat

Angle = rad2deg(atan(1.5/4096));
PIVW_rot = imrotate(IM,Angle);

%% Second: affine transformation to rectify image
% Correct PIV Surf Camera Angle
U1 = [392 1179 ; 379 3000 ; 4045 3002 ; 4029 1180]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [386 1179.5 ; 386 3001; 4037 3001; 4037 1179.5]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVWater_CamAngle =imwarp(PIVW_rot,T1,'cubic');

%% Third: possible other rototranslation to match the images
