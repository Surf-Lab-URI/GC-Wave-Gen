function [PIVAir_CamAngle] = PIVAir_CamAngle_Correction(IM)

%% PIV Air CamAngle Correction

%% First: we checked the surface is flat

%% Second: affine transformation to rectify image
% Correct PIV Camera Angle
U1 = [384 2161 ; 361 57 ; 4017.5 57 ; 3997 2161]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [ (384+361)/2 2161 ; (384+361)/2 57 ; (4017.5+3997)/2 57 ; (4017.5+3997)/2 2161 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
PIVAir_CamAngle =imwarp(IM,T1,'cubic');

%% Third: possible other rototranslation to match the images
