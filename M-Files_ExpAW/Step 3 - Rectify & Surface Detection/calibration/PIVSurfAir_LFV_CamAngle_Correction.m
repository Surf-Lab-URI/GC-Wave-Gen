function [PIVSurfA_CamAngle] = PIVSurfAir_LFV_CamAngle_Correction(PIVSurfA_Undistorted)

%% PIVSurf Air - LFV CamAngle Correction

%% First: rotation to make surface flat
Angle = rad2deg(atan((3166.5-3164.5)/(7865-51)));
PIVSurfA_rot = imrotate(PIVSurfA_Undistorted,Angle);

%% Second: affine transformation to rectify image
% Correct PIV Surf Camera Angle
% U1 = [ 7217 2402 ; 7351 450 ; 66 372 ; 170 2395]; % The coordinates of 
U1 = [2536 3120 ; 2477 60 ; 4955 60 ; 4907 3120]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
% X1 =  [ 7351 2395 ; 7351 450 ; 121 450 ; 121 2395]; % The coordinates of 
X1 = [ (2536+2477)/2 3120 ; (2536+2477)/2 60 ; (4955 + 4907)/2 60 ; (4955 + 4907)/2 3120 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVSurfA_CamAngle =imwarp(PIVSurfA_rot,T1,'cubic');

%% Third: possible other rototranslation to match the images