function [scaled,scaledSmallCrop,scaledCroppedToPIV] = SurfImgToPIVDims(imgPivsurf)
%Surface detection and Creating Masks
U1 = [147 49;2024 57; 1995 1004; 161 999];
X1 = [147 49; 2024 49; 2024 1004; 147 1004];
T1 = fitgeotrans(U1,X1,'projective');
scaled=imwarp(imgPivsurf,T1,'cubic');

scaled=imresize(scaled,176.9769/105.5880);%Resizing to match PIV resolution but not dimensions

scaledSmallCrop = scaled(30:3525, 100:end-200);%cropping to eliminate artifacts of lens distortion correction around the edges.

scaledCroppedToPIV=scaled(30:3525,755:755+2047); %cropping to match PIV image dimensions.

end

