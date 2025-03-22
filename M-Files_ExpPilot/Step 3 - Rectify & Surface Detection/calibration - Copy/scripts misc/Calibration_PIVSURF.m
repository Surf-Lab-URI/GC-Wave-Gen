% List of commands used for PIVSURF calibration

clear
close all
clc

%% Combo grid
% Construct Combo_grid (image composed by grid, flat surface and
% plumblines)
IMgrid = fliplr(load_Image_IOCoreView_12MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\PIV_match\05232022_plumbline_Movie1_Scene3_PIVSURF_0.raw'));
IMflat = fliplr(load_Image_IOCoreView_12MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\flat\Movie1_Scene1_PIVSURF_00.raw'));
IMflat(IMflat>20)=10;
IMpl = fliplr(load_Image_IOCoreView_12MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\plumbline\Movie1_Scene1_PIVSURF_0.raw'));
Combo_grid = IMgrid;
Combo_grid(1501:1520,:) = 5*double(imadjust(uint8(IMflat(1501:1520,:))));
Combo_grid(1:1000,121:210) = IMpl(1:1000,121:210);
Combo_grid(1:1000,3915:4000) = IMpl(1:1000,3915:4000);

% Compute Combo_grid undistorted
[Combo_grid_lens] = CorrectPIVSURFLensDistortion(Combo_grid);
figure;imagesc(Combo_grid_lens);colormap gray

% Compute Combo_grid with camera angle correction to check straight lines,
% vertical plumblines and horizontal flat surface
[Combo_grid_CamAngle_Corrected] = Correct_Angle_PIVSURF(Combo_grid_lens);
figure;imagesc(Combo_grid_CamAngle_Corrected)

%% IMgrid used from PIVSURF to match PIV resolution
%%% Without polynomial
% % % IMgrid_lens = CorrectPIVSURFLensDistortion(IMgrid);
% % % IMgrid_CamAngle_Corrected = fliplr(Correct_Angle_PIVSURF(IMgrid_lens));
%%% With polynomial
PIVSURF_grid_rect = load('rect.mat','PIVSURF_grid_rect');
PIVSURF_grid_rect = fliplr(PIVSURF_grid_rect.PIVSURF_grid_rect);

% Find checkerboard corners
[imagePoints, boardSize] = detectCheckerboardPoints(uint16(PIVSURF_grid_rect));
imagePoints(1:9*(boardSize(1)-1),:) = [];
PIVSURF_iP = [];
% Reshape to match PIV points
for ii = 1:boardSize(1)-1
    for i = 1:boardSize(2)-1-9
        PIVSURF_iP((boardSize(2)-1-9)*(ii-1)+i,:) = imagePoints(((boardSize(1)-1)+1-ii)+(boardSize(1)-1)*(i-1),:);
    end
end
% Delete exceeding points
PIVSURF_iP(1:boardSize(2)-1-9,:) = [];
figure;imagesc(PIVSURF_grid_rect)
hold on;plot(PIVSURF_iP(:,1),PIVSURF_iP(:,2),'ro');

% Points used: 1) O-9 (O left-under) ; 2) B-9 (B left-under) ; 3) B-22 (B left-up) ; 4) O-22 (O left-up)
Up = [imagePoints(cursor_info2(1).DataIndex,:) ; imagePoints(cursor_info(1).DataIndex,:);...
      imagePoints(cursor_info(2).DataIndex,:) ; imagePoints(cursor_info2(2).DataIndex,:)];
  
% Save matching points
PIVSURF_iP2 = PIVSURF_iP; % this is to save points on polynomial image IMgrid_CamAngle_Corrected.IMgrid_CamAngle_Corrected
save('\\spray1-10g\d\Shoaling_waves\M-files\Step 3 - Rectify & Surface Detection\calibration\PIV_PIVsurf_matching_points', '-append', 'Up', 'PIVSURFgrid_CamAngle_Corrected', 'PIVSURF_iP', 'PIVSURF_iP2')