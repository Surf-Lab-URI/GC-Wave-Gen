function [PIVWater_CamAngle] = PIVWater_CamAngle_Correction(IM)

%% PIV Air CamAngle Correction

%% First: no need of further rotation to make surface flat

%% Second: affine transformation to rectify image
% Correct PIV Surf Camera Angle
U1 = [428 720 ; 403 3005 ; 3965.5 3005 ; 3933.75 720]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [ (428+403)/2 720 ; (428+403)/2 3005 ; (3965.5+3933.75)/2 3005 ; (3965.5+3933.75)/2 720 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVWater_CamAngle =imwarp(IM,T1,'cubic');

%% Third: possible other rototranslation to match the images
