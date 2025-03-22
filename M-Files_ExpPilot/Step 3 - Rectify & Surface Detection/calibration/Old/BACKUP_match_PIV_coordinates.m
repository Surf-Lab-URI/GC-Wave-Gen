%% Transform PIV Water measurements in PIVSurf Water coordinates

%% Step 1: Inverse transform of the SURFACE - From PIV Water to PIVSurf Water
% This step is necessary to verify the matching between original surface and the transformed one 

%%% Transformation from PIV Water to PIVSurf Water coordinate system of the
%%% surface
% Points to use for projection of Landmarks on the PIVSurf Water grid
Rp1 = fliplr([ 327 3058 ; 731 3056 ; 1135 3052 ; 1541 3048 ; 1944 3045 ]); % GRID POINTS on the first straight line RIGHT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_R = -3.25; DY_R = 404.25; % main DIFFERENCE in the horizontal (DX_R) and in the vertical (DY_R) 
Lp1 = [1013 311 ; 1010 714 ; 1007 1118 ; 1004 1523 ; 1000 1927 ]; % GRID POINTS on the first straight line LEFT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_L = -3.25; DY_L = 404; % main DIFFERENCE in the horizontal (DX_L) and in the vertical (DY_L) 
% Landmarks in the PIVSurf Water grid 
Xp1(2:3,1:2) = [Rp1(end,1)+DX_R Rp1(end,2)+DY_R ; Rp1(end,1)+4*DX_R Rp1(end,2)+4*DY_R];
Xp1([1,4],1:2) = [Lp1(end,1)+DX_L Lp1(end,2)+DY_L ; Lp1(end,1)+4*DX_L Lp1(end,2)+4*DY_L];

U1 = Xp1; % COORDINATES IN PIVSURF WATER
X1 = [ 404 1195 ; 3447 1198 ; 3447 3018 ; 405 3019 ]; % COORDINATES IN PIV WATER

T1 = maketform('projective',U1,X1); %Transformation from PIVSurf Water to PIV Water

% Eliminate the rototranslation added after the transformation
Yinv2 = PIVW_PIVSurfW1_Surface-DY;
Yinv2 = Minv*[XPIVW_PIVSurfW1_Surface;Yinv2];
Xinv = XPIVW_PIVSurfW1_Surface;
Yinv = Yinv2(2,:);

% Inverse mapping of the surface
[Uinv2,Vinv2] = tforminv(T1,Xinv,Yinv);
Uinv = XPIVSurfW1_Surface;
Vinv = interp1(Uinv2,Vinv2,Uinv,'linear','extrap');

%%% Original surface: XPIVSurfW1_Surface,PIVSurfW1_Surface
%%% Transformed surface: Uinv, Vinv
% figure;imagesc(PIVSurfW1_CamAngle);colormap gray;caxis([0 200]);hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')
% figure;plot(Uinv,PIVSurfW1_Surface-Vinv,'r')

%% Step 2: Inverse transform of the PIV Water VELOCITY - From PIV Water to PIVSurf Water

%%% Inverse transformation (from PIV Water to PIVSurf Water)
% Points to use for projection of Landmarks on the PIVSurf Water grid
Rp1 = fliplr([ 327 3058 ; 731 3056 ; 1135 3052 ; 1541 3048 ; 1944 3045 ]); % GRID POINTS on the first straight line RIGHT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_R = -3.25; DY_R = 404.25; % main DIFFERENCE in the horizontal (DX_R) and in the vertical (DY_R) 
Lp1 = [1013 311 ; 1010 714 ; 1007 1118 ; 1004 1523 ; 1000 1927 ]; % GRID POINTS on the first straight line LEFT of the tape (calibration image calibration_ExpPilot_Scene4_PIVSURF Water_11.raw) 
DX_L = -3.25; DY_L = 404; % main DIFFERENCE in the horizontal (DX_L) and in the vertical (DY_L) 
% Landmarks in the PIVSurf Water grid 
Xp1(2:3,1:2) = [Rp1(end,1)+DX_R Rp1(end,2)+DY_R ; Rp1(end,1)+4*DX_R Rp1(end,2)+4*DY_R];
Xp1([1,4],1:2) = [Lp1(end,1)+DX_L Lp1(end,2)+DY_L ; Lp1(end,1)+4*DX_L Lp1(end,2)+4*DY_L];

U1 = Xp1; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X1 = [ 404 1195 ; 3447 1198 ; 3447 3018 ; 405 3019 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.

T2inv = maketform('projective',X1,U1);

%%% Interpolate PIVWater measurements with pixel resolution
[XX,YY] = meshgrid(4:4:size(CompVelWater.INTdelx,2)*4,4:4:size(CompVelWater.INTdelx,1)*4);
[XXq,YYq] = meshgrid(1:size(PIV1_W,2),1:size(PIV1_W,1));
% VVq = griddata(XX,YY,CompVelWater3.INTdelx,XXq,YYq);
VVq = interp2(XX,YY,CompVelWater.INTdelx,XXq,YYq);

%%% Rototranslate the image
DY = 17; % Same translation used in Extract_PIVSurf Water
RotAngle = -39/6107; % Same rotation angle used in Extract_PIVSurf Water
Minv = [cos(-RotAngle) sin(-RotAngle); -sin(-RotAngle) cos(-RotAngle)];

VVq = VVq; %.*Mask1_W;
% VVVq = VVq; VVVq(isnan(VVq)) = 9999;
IMrot = imrotate(VVq,-rad2deg(RotAngle));
VVVq = IMrot; VVVq(isnan(IMrot)) = 9999;
IMtr = imtranslate(VVVq,[0 -45]);
IMtr(IMtr==9999) = NaN;
%IMrot = imrotate(IMtr,-rad2deg(RotAngle));

% [CompVelW_inv,Xpos,Ypos] = imtransform(IMrot,T2inv,'XYScale',1);
[CompVelWater_inv,XposInv,YposInv] = imtransform(IMtr,T2inv,'XYScale',1);

% Water velocity with pixel resolution in PIVSurf Water coordinates
[XX2,YY2] = meshgrid(XposInv(1):XposInv(end),YposInv(1):YposInv(end));
[XXq2,YYq2] = meshgrid(1:XposInv(end),1:YposInv(end));
% CompVelWater_PIVSurfW = griddata(XX2,YY2,CompVelWater_inv,XXq2,YYq2);
CompVelWater_PIVSurfW = interp2(XX2,YY2,CompVelWater_inv,XXq2,YYq2);

%%% Transformed PIV Water velocity: CompVelW_inv
% figure;imagesc(XposInv(1):XposInv(end),YposInv(1):YposInv(end),CompVelW_inv);
% hold on;plot(XPIVSurfW1_Surface,PIVSurfW1_Surface,'r')
% plot(Uinv,Vinv,'g')

%% Step 3: Transform surface and PIV Water velocity - from PIVSurf Water to PIV Air coordinates
U2 = [999 1927 ; 1009 816 ; 3056 832 ; 3045 1944 ]; % The coordinate
% of four points in the PIVSurf Air - LFV surface image (the grid calibration image).
X2 = [481 1877 ; 490 44 ; 3825 63 ; 3815 1894 ]; % The coordinates of  the four
% points in the PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T3 = maketform('projective',U2,X2);

[Xsurf2,Ysurf2] = tformfwd(T3,Uinv,Vinv);
Xsurf = round(Xsurf2(1):Xsurf2(end));
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf);

%%% This is the actual resized image; we can retrieve it for checking
[CompVelWater_PIVA2,Xpos,Ypos] = imtransform(CompVelWater_PIVSurfW,T3,'XYScale',1);
CompVelWater_PIVA2(:,[1:1220,5732:end]) = NaN;
CompVelWater_PIVA2([1:2560,5889:end],:) = NaN;

% Water velocity with pixel resolution in PIVSurf Water coordinates
[XX2,YY2] = meshgrid(Xpos(1):Xpos(end),Ypos(1):Ypos(end));
[XXq2,YYq2] = meshgrid(1:Xpos(end)-1,1:Ypos(end));
CompVelWater_PIVA = interp2(XX2,YY2,CompVelWater_PIVA2,XXq2,YYq2);

%% Warp CompVelWater to match PIV_Surface (PIVSurf Air in PIV Air coordinates)
movingPoints = Ysurf(1103:5710);
fixedPoints = PIV_LFV_Surface(4333:8940);
DeltaY = round(fixedPoints-movingPoints);
Xmov = 1:size(CompVelWater_PIVA,1);
CompVelWater_PIVA_warp = nan(size(PIV1_A,1),size(PIV1_A,2));
for i = 1:length(DeltaY)
    Xintrp = 1:size(CompVelWater_PIVA,1)/(size(CompVelWater_PIVA,1)+DeltaY(i)):size(CompVelWater_PIVA,1);
    Dumb = interp1(Xmov,CompVelWater_PIVA(:,i),Xintrp);
    if DeltaY(i)>0
        Xd = DeltaY(i):length(Dumb);
        Xc = 1:4700;
    elseif DeltaY(i)==0
        Xd = DeltaY(i)+1:length(Dumb);
        Xc = 1:4700;
    else
        Xd = 1:length(Dumb);
        Xc = 1:4700+DeltaY(i);
    end
    CompVelWater_PIVA_warp(Xc,i) = Dumb(Xd:Xd+min(length(Dumb),4700)-1);
end

[PIV_Mask_W1] = PIVWater_Mask(CompVelWater_PIVA, PIV_LFV_Surface(4333:8940));
CompVelWater_PIVA_warp = CompVelWater_PIVA.*PIV_Mask_W1;

CompVelWater_PIVAir.delta_x = CompVelWater_PIVA_warp(4:4:end-4);