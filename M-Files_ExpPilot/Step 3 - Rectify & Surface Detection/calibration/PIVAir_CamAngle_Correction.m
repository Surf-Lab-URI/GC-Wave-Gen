function [PIVAir_CamAngle] = PIVAir_CamAngle_Correction(IM)

%% PIV Air CamAngle Correction

%% First: rotation to make surface flat

Angle = rad2deg(atan(4/4096));
PIVA_rot = imrotate(IM,-Angle);

%% Second: affine transformation to rectify image
% Correct PIV Camera Angle
U1 = [473 2140 ; 455 59 ; 3807 55 ; 3792 2140]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [464 2140 ; 464 57; 3799.5 57; 3799.5 2140]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVAir_CamAngle =imwarp(PIVA_rot,T1,'cubic');

%% Third: possible other rototranslation to match the images
