%% Backup on calibration steps

%% 1) PIV Air
%% PIV flat surface and plumblines
% Flat surface
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\flat_surface\flat_surface_Scene9\RAW\PIV Air\';
filename = [LoadPath 'flat_surface_Scene9_PIV Air_009.raw'];
[IMflat] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(IMflat);colormap gray;caxis([0 100])

% Plumblines
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\plumblines\plumblines_Scene4\RAW\PIV Air\';
filename = [LoadPath 'plumblines_Scene4_PIV Air_22.raw'];
[IMplumb] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(IMplumb);colormap gray;caxis([0 100])

% Combo image with flat surface and plumblines
IMcombo = IMflat;
IMcombo(:,[301:450,2031:2180,3951:end]) = IMplumb(:,[301:450,2031:2180,3951:end]);
% figure;imagesc(IMcombo);colormap gray;caxis([0 100])

% Points used for PIV Air CamAngle correction
U1 = [384 2161 ; 361 57 ; 4017.5 57 ; 3997 2161];
X1 = [ (384+361)/2 2161 ; (384+361)/2 57 ; (4017.5+3997)/2 57 ; (4017.5+3997)/2 2161 ];
T1 = fitgeotrans(U1,X1,'projective'); 
IMcombo_CamAngle = imwarp(IMcombo,T1,'cubic');

%% Calibration point to refer PIV coordinates
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\calibration_grid\calibration_grid\calibration_grid_Scene4\RAW\PIV Air\';
filename = [LoadPath 'calibration_grid_Scene4_PIV Air_20.raw'];
[IM1] = (load_Image_IOCoreView_12MP(filename));
IM1_CamAngle = imwarp(IM1,T1,'cubic');
figure;imagesc(IM1_CamAngle);colormap gray;clim([0 100])

% Points used for PIV coordinates references
X2 = [ 595 2033 ; 616 296 ; 3781 322; 3768 2062 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------------------------------------- % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2) PIVSurf Air - LFV
%% Undistort image
% % % use CameraCalibrator here to find intrinsic camera parameters % % %
load cameraParams.mat

%% PIVSurf Air - LFV flat surface and plumblines
% Flat surface
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\flat_surface\flat_surface_Scene4\RAW\PIVSurf Air - LFV\';
imagename = [LoadPath 'flat_surface_Scene4_PIVSurf Air - LFV_07.raw'];
[LFVflat] = load_Image_IOCoreView_48MP(imagename);
LFVflat_Undistorted = undistortImage(LFVflat,cameraParams);
% figure;imagesc(LFVflat_Undistorted);colormap gray;clim([0 300])

% Plumblines
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\plumblines\plumblines_Scene4\RAW\PIVSurf Air - LFV\';
imagename = [LoadPath 'plumblines_Scene4_PIVSurf Air - LFV_35.raw'];
[LFVplumb] = load_Image_IOCoreView_48MP(imagename);
LFVplumb_Undistorted = undistortImage(LFVplumb,cameraParams);
% figure;imagesc(LFVplumb_Undistorted);colormap gray;clim([0 300])

% Combo surface
LFVcombo = LFVflat_Undistorted;
LFVcombo(:,[2440:2560,3600:3670,4880:4980]) = LFVplumb_Undistorted(:,[2440:2560,3600:3670,4880:4980]);
% figure;imagesc(LFVcombo);colormap gray;clim([0 300])

% Rotation to make surface flat
Angle = rad2deg(atan((3166.5-3164.5)/(7865-51)));
LFVcombo_rot = imrotate(LFVcombo,Angle);
% figure;imagesc(LFVcombo_rot);colormap gray;clim([0 300])

% Points used for PIVSurf Air - LFV CamAngle correction (flat surface and
% straight plumblines)
U1 = [2536 3120 ; 2477 60 ; 4955 60 ; 4907 3120];
X1 = [ (2536+2477)/2 3120 ; (2536+2477)/2 60 ; (4955 + 4907)/2 60 ; (4955 + 4907)/2 3120 ];
T1_LFV = fitgeotrans(U1,X1,'projective'); 
LFVcombo_CamAngle = imwarp(LFVcombo_rot,T1_LFV,'cubic');
figure;imagesc(LFVcombo_CamAngle);colormap gray;clim([0 300])

%% Calibration point to refer to PIV Air coordinates
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\calibration_grid\calibration_grid\calibration_grid_Scene1\RAW\PIVSurf Air - LFV\';
imagename = [LoadPath 'calibration_grid_Scene1_PIVSurf Air - LFV_20.raw'];
[LFVgrid] = load_Image_IOCoreView_48MP(imagename);
LFVgrid = undistortImage(LFVgrid,cameraParams);
LFVgrid = imrotate(LFVgrid,Angle);
LFVgrid =  imwarp(LFVgrid,T1_LFV,'cubic');
figure;imagesc(LFVgrid);colormap gray;clim([0 300])

% Points used for PIV coordinates references
U2 = [ 2872 3085 ; 2883 1890 ; 5051 1913 ; 5038 3105];
T2 = maketform('projective',U2,X2);
[Resized_LFVgrid,XPos,YPos] =  imtransform(LFVgrid,T2,'XYScale',1);
figure;imagesc(XPos(1):XPos(2),YPos(1):YPos(2),Resized_LFVgrid);colormap gray;axis([1,4143,1,3088]);
hold on;plot(X2(:,1),X2(:,2),'ro')

% NOTE!!! The previous code does not conserve plumbline verticality; if
% surface not matching well with PIVAir, try with the following rows (imresize and imtranslate)
% Scale = nanmean(nanmean(dist(X2')./dist(U2')));
% LFVRes = imresize(LFVgrid,Scale);
% % [3768 2062] in PIVAir corresponds to [7356 4533] in LFVRes = imresize(LFVgrid,Scale). Thus LFVRes = imtranslate(LFVRes,-[7356-3768,4533-2062]);
% LFVRes = imtranslate(LFVRes,-[7356-3768,4533-2062]);
% figure;imagesc(LFVRes);colormap gray;axis([1,4143,1,3088])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------------------------------------- % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3) PIVWater 
%% PIVWater flat surface and plumblines
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\flat_surface\flat_surface_Scene9\RAW\PIV Water\';
filename = [LoadPath 'flat_surface_Scene9_PIV Water_009.raw'];
[IMflat] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(IMflat);colormap gray;clim([0 100])

% Plumblines
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\plumblines\plumblines_Scene4\RAW\PIV Water\';
filename = [LoadPath 'plumblines_Scene4_PIV Water_22.raw'];
[IMplumb] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(IMplumb);colormap gray;clim([0 100])

% Combo image with flat surface and plumblines
IMcombo = IMflat;
IMcombo(:,[341:500,2031:2100,3880:3990]) = IMplumb(:,[341:500,2031:2100,3880:3990]);
% figure;imagesc(IMcombo);colormap gray;clim([0 100])

% Points used for PIV Air CamAngle correction
U1 = [428 720 ; 403 3005 ; 3965.5 3005 ; 3933.75 720];
X1 = [ (428+403)/2 720 ; (428+403)/2 3005 ; (3965.5+3933.75)/2 3005 ; (3965.5+3933.75)/2 720 ];
T1 = fitgeotrans(U1,X1,'projective'); 
IMcombo_CamAngle = imwarp(IMcombo,T1,'cubic');
figure;imagesc(IMcombo_CamAngle);colormap gray;clim([0 200])

%% Calibration point to refer PIV coordinates
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\calibration_grid\calibration_grid\calibration_grid_Scene4\RAW\PIV Water\';
filename = [LoadPath 'calibration_grid_Scene4_PIV Water_20.raw'];
[IM1] = (load_Image_IOCoreView_12MP(filename));
IM1_CamAngle = imwarp(IM1,T1,'cubic');
figure;imagesc(IM1_CamAngle);colormap gray;clim([0 100])

% Points used for PIV Water coordinates references
X2_W = [ 630 875 ; 610 2715 ; 3706 2743; 3718 903 ];
% hold on;plot(X2_W(:,1),X2_W(:,2),'ro')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------------------------------------- % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4) PIVSurf Water

%% PIVSurf Water flat surface and plumblines
% Flat surface
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\flat_surface\flat_surface_Scene9\RAW\PIVSurf Water\';
filename = [LoadPath 'flat_surface_Scene9_PIVSurf Water_010.raw'];
[PIVSurfWflat] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(PIVSurfWflat);colormap gray;caxis([0 100])

% Plumblines
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\plumblines_and_flat_surface\plumblines\plumblines_Scene4\RAW\PIVSurf Water\';
filename = [LoadPath 'plumblines_Scene4_PIVSurf Water_22.raw'];
[PIVSurfWplumb] = (load_Image_IOCoreView_12MP(filename));
% figure;imagesc(PIVSurfWplumb);colormap gray;caxis([0 100])

% Combo image with flat surface and plumblines
PIVSurfWcombo = PIVSurfWflat;
PIVSurfWcombo(:,[1011:1070,2065:2110,3221:3320]) = PIVSurfWplumb(:,[1011:1070,2065:2110,3221:3320]);
% figure;imagesc(PIVSurfWcombo);colormap gray;caxis([0 100])

% Points used for PIV Air CamAngle correction
U1_W = [1052 2010 ; 1040 951 ; 3280.5 951 ; 3265 2010];
X1_W = [ (1052+1040)/2 2010 ; (1052+1040)/2 951 ; (3280.5+3265)/2 951 ; (3280.5+3265)/2 2010 ];
T1_PIVsurfW = fitgeotrans(U1_W,X1_W,'projective'); 
PIVSurfWcombo_CamAngle = imwarp(PIVSurfWcombo,T1_PIVsurfW,'cubic');
% figure;imagesc(PIVSurfWcombo_CamAngle);colormap gray;caxis([0 100])

%% Calibration point to refer PIV Water coordinates
LoadPath = '\\spray3\d\data\EXPERIMENTS\Calibration\calibration_grid\calibration_grid\calibration_grid_Scene4\RAW\PIVSURF Water\';
filename = [LoadPath 'calibration_grid_Scene4_PIVSurf Water_20.raw'];
[PIVSurfWgrid] = (load_Image_IOCoreView_12MP(filename));
PIVSurfWgrid_CamAngle = imwarp(PIVSurfWgrid,T1_PIVsurfW,'cubic');
figure;imagesc(PIVSurfWgrid_CamAngle);colormap gray;clim([0 100])

% Points used for PIV Water coordinates references
%%% Projection of the grid points under the water
% Left points LP1 and LP2
LP_Inch8th = [1203 1963 ; 1203 1915 ; 1204 1866 ; 1204 1819 ; 1204 1772 ; 1205 1724 ; 1205 1675 ; 1206 1627 ; 1206 1579 ; 1207 1531 ; 1207 1483 ; 1207 1436 ; 1208 1387 ; 1208 1338 ; 1208 1290 ];
DeltaY_LP = -mean(diff(LP_Inch8th(:,2)));
CC1 = polyfit(LP_Inch8th(:,2),LP_Inch8th(:,1),1); % straight line fitting landmarks
LP1(1,2) = LP_Inch8th(1,2)+4*DeltaY_LP;
LP1(1,1) = polyval(CC1,LP1(1,2));
LP2(1,2) = LP1(1,2)+(3*8)*DeltaY_LP;
LP2(1,1) = polyval(CC1,LP2(1,2));
hold on;plot(LP_Inch8th(:,1),LP_Inch8th(:,2),'cx',polyval(CC1, 1:6000),1:6000,'g');plot(LP1(1),LP1(2),'ro',LP2(1),LP2(2),'rx')
% Right points RP1 and RP2
RP_Inch8th = [ 3163 1986 ; 3163 1938 ; 3164 1890 ; 3164 1841 ; 3164 1793 ; 3165 1745 ; 3165 1697 ; 3165 1649 ; 3166 1601 ; 3166 1553 ; 3166 1503 ; 3167 1455 ;  3167 1407 ; 3168 1359 ;  3168 1311];
DeltaY_RP = -mean(diff(RP_Inch8th(:,2)));
CC2 = polyfit(RP_Inch8th(:,2),RP_Inch8th(:,1),1); % straight line fitting landmarks
RP1(1,2) = RP_Inch8th(1,2)+4*DeltaY_RP;
RP1(1,1) = polyval(CC2,RP1(1,2));
RP2(1,2) = RP1(1,2)+(3*8)*DeltaY_RP;
RP2(1,1) = polyval(CC2,RP2(1,2));
hold on;plot(RP_Inch8th(:,1),RP_Inch8th(:,2),'cx',polyval(CC2, 1:6000),1:6000,'g');plot(RP1(1),RP1(2),'ro',RP2(1),RP2(2),'rx')

U2_W = [ LP1 ; LP2 ; RP2 ; RP1];
T2_W = maketform('projective',U2_W,X2_W);
[Resized_PIVSurfWgrid_CamAngle,XPos_W,YPos_W] =  imtransform(PIVSurfWgrid_CamAngle,T2_W,'XYScale',1);
figure;imagesc(XPos_W(1):XPos_W(2),YPos_W(1):YPos_W(2),Resized_PIVSurfWgrid_CamAngle);colormap gray;axis([1,4143,1,3088]);
hold on;plot(X2_W(:,1),X2_W(:,2),'ro')

%% Calibration points for PIV coordinates
% Transformation from PIVSurf Water to PIV Air
U2_W2 = [ 1203 1963 ; 1213 1001 ; 3171 1016; 3163 1986];
X2_W2 = [ 595 2033 ; 615 454 ; 3780 480; 3768 2062 ];
T2_W2 = maketform('projective',U2,X2);

%% 5) From PIVWater to PIVAir
%1) PIVWater --> inverse T2_W --> PIVSurf Water (both surface and velocity)
%2) PIVSurf Water --> T2_W2 --> PIV Air
%3) Then remember to warp velocity to match PIVSurf Air
