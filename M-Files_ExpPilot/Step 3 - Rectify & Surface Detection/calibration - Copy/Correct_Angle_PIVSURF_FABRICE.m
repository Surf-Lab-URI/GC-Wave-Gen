function [PIVSurf_Corrected] = Correct_Angle_PIVSURF_FABRICE (PIVSurf_CR1)
% The CorrectPIVSurf function corrects the camera angle first, and then, it
% resizes the image to have the same resolution of the fused PIV images. It
% is important to notice that the PIV surface images need  to have the same
% resolution of the fused PIV image. In order to obtain the same resolution
% as fused PIV image, the "imtransform" function is used.

% %Rotate a bit for flat surface to be horizontal
% PIVSurf_CR1=imrotate(PIVSurf_CR1,-asin(4/1792)/pi*180,'bicubic');

% % % U1 = []; % grid for matching PIV coordinate
% % % X1 = fliplr([5812 6971; 120 6971; 120 274; 5812 274]); % from PIV
% coordiubate

% Correct PIV Surf Camera Angle
U1 = [11 200; 49 914; 1995 914; 2034 200]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [11 200; 11 914; 2034 914; 2034 200]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.
T1 = fitgeotrans(U1,X1,'projective'); 
CamAngle =imwarp(PIVSurf_CR1,T1,'cubic');

% Resize PIV Surf images to the Size of Fused PIV images
%Tweek = 7;
U2 = [473 655; 474 941; 1784 939; 1783 652]; % The coordinate
% of four points in the PIV surface image (the grid calibration image).
X2 = [1356 1003; 1363 2179; 6655 2167; 6647 997]; % The coordinates of  the four
% points in the fused PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T2 = maketform('projective',U2,X2);
[Resized_PIVSurf,XPos,YPos] =  imtransform(CamAngle,T2,'XYScale',1);

%one more rotation after the transformation
Resized_PIVSurf=imrotate(Resized_PIVSurf,-asin(1/1792)/pi*180,'bicubic');

PIVSurf_Corrected.img  = Resized_PIVSurf;
PIVSurf_Corrected.Xpos = XPos;
PIVSurf_Corrected.Ypos = YPos;

% XPos and YPos are a two-element, real vector that  together  specifiy the
% spatial location of the output image B in the 2D output space XY. The two
% elements of XPos and YPos give the x-coordinates and y-coordinates of the
% first and last columns of B, respectively.

%Resized_PIVSurf(round(abs(YPos(1))):round(abs(YPos(1)))+s(1)-1,round(abs(XPos(1))):round(abs(XPos(1)))+s(2)-1)=FusedPIV;
end