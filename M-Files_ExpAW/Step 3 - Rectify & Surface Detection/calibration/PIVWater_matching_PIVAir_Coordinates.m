function [PIVWater_CRR] = PIVWater_matching_PIVAir_Coordinates(PIVWater_CamAngle)

%% Find straight line which best approximates the position of landmarks in PIVSurf Water with PIV coordinates
%%% Points [Xp1 Yp1] and [Xp2 Yp2] are retrieved from Resized_PIVSurfW'
% Points on the left side of the image
Xp1 = [ -788.1 -124.1 542.8 1210.9 1875.9 ];
Yp1 = [ 489.4 485.4 484.4 483.4  478.4 ];
A1 = polyfit(Xp1,Yp1,1);
% Points on the right side of the image
Xp2 = [ -775.1 -109.1 558 1226.9 1892.9 ];
Yp2 = [ 3824.4 3821.4 3819.4  3816.4 3812.4 ];
A2 = polyfit(Xp2,Yp2,1);

%% Find landmarks of PIV Water in transformed coordinates
% Point 1 (upper left)
X2(1,1) = Xp1(end)+1*mean(diff(Xp1));
X2(1,2) = polyval(A1,X2(1,1));
% Point 2 (lower left)
X2(2,1) = Xp1(end)+4*mean(diff(Xp1));
X2(2,2) = polyval(A1,X2(2,1));
% Point 3 (lower right)
X2(3,1) = Xp2(end)+4*mean(diff(Xp2));
X2(3,2) = polyval(A2,X2(3,1));
% Point 4 (upper right)
X2(4,1) = Xp2(end)+1*mean(diff(Xp2));
X2(4,2) = polyval(A2,X2(4,1));

%%% Check if it is right
% figure;imagesc(YPosW(1):YPosW(2),XPosW(1):XPosW(2),Resized_PIVSurfW');colormap gray;caxis([0 50])
% ylim([1 4140]);xlim([0 3090])
% hold on;plot(-1000:5000,polyval(A1,-1000:5000),'r',-1000:5000,polyval(A2,-1000:5000),'m')
% plot(X2(1,1),X2(1,2),'rx',X2(2,1),X2(2,2),'bx',X2(3,1),X2(3,2),'gx',X2(4,1),X2(4,2),'cx')

%% Transformation to match PIV Water and PIV Air
% Correct PIV Surf Camera Angle
U1 = [404 1195 ; 404 3019 ; 3447 3019 ; 3447 1198 ]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = fliplr(X2); % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.

T2 = maketform('projective',U1,X1);
[Resized_PIVWater,XPos,YPos] =  imtransform(PIVWater_CamAngle,T2,'XYScale',1);

PIVWater_CRR.img = Resized_PIVWater;
PIVWater_CRR.Xpos = XPos;
PIVWater_CRR.Ypos = YPos;
PIVWater_CRR.Tform = T2;
