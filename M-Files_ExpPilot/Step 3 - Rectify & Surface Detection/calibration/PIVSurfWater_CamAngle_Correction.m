function [PIVSurfW_CamAngle] = PIVSurfWater_CamAngle_Correction(IM)

%% PIVSurf Air - LFV CamAngle Correction

%% First: rotation to make surface flat
Angle = rad2deg(atan(6/4096));
PIVSurfW_rot = imrotate(IM,-Angle);

%% Second: affine transformation to rectify image
% Correct PIV Surf Camera Angle
U1 = [ 172 1943 ; 3821 1953 ; 3849 338 ; 145 329]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 =  [ 158.5 1948 ; 3835 1948 ; 3835 333.5 ; 158.5 333.5]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVSurfW_CamAngle =imwarp(PIVSurfW_rot,T1,'cubic');

%% Third: possible other rototranslation to match the images