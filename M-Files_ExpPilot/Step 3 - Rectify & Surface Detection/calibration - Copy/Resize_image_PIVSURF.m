function [PIVSurf_Corrected,T2_bis] = Resize_image_PIVSURF(PIV2_CamAngle_Corrected)

%% Resize PIV Surf images to the Size of PIV images
load PIV_match.mat U1 X1 Xt
T1 = fitgeotrans(U1,X1,'polynomial',4);
PIV2_CamAngle_Corrected(isnan(PIV2_CamAngle_Corrected)) = 0;
PIV2_CamAngle_Corrected = imwarp(PIV2_CamAngle_Corrected,T1,'linear');
PIV2_CamAngle_Corrected = imtranslate(PIV2_CamAngle_Corrected, [(X1(1,1)-Xt(1,1)) (X1(1,2)-Xt(1,2))] );
PIV2_CamAngle_Corrected = fliplr(PIV2_CamAngle_Corrected); 

% Load matching points
P = load('PIV_PIVsurf_matching_points.mat','Xp');
% % % Up = round(P.Up); % The coordinate
Up2 = [-X1(:,1)+size(PIV2_CamAngle_Corrected,2)+1,X1(:,2)];
Up2([1:26,end-8*26+1:end],:) = [];
Up2 = flipud(reshape(permute(reshape(Up2,[26,31,2]),[2,1,3]),26*31,2));
Up = [ Up2(end-30,:) ; Up2(end,:) ; Up2(31,:) ; Up2(1,:) ]; % The coordinate
% of four points in the PIV surface image (the grid calibration image).
Xp = P.Xp; % The coordinates of  the four
% points in the fused PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.
T2_bis = maketform('projective',Up,Xp);
[Resized_PIVSurf_bis,XPos,YPos] =  imtransform(PIV2_CamAngle_Corrected,T2_bis,'XYScale',1);

%%% I use a polynomial transformation to map all the matching points
%%% between PIVSURF and PIV. Then, I find the origin of the Resized_PIVSurf
%%% (Pos) and from that I can retrieve XPos and YPos

Resized_PIVSurf = Resized_PIVSurf_bis;

PIVSurf_Corrected.img  = Resized_PIVSurf;
PIVSurf_Corrected.XPos = XPos;
PIVSurf_Corrected.YPos = YPos;

% XPos and YPos are a two-element, real vector that  together  specifiy the
% spatial location of the output image B in the 2D output space XY. The two
% elements of XPos and YPos give the x-coordinates and y-coordinates of the
% first and last columns of B, respectively.