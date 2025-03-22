% List of commands used for PIV calibration

clear
close all
clc

%% Combo grid
% Construct Combo_grid (image composed by grid, flat surface and
% plumblines)
IMgrid = fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\PIVSURF_match\Movie1_Scene1_PIV_00.raw'));
IMflatAfter2 = fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\flat\05232022_flat_Movie1_Scene1_PIV_000.raw'));
IMpl = fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\plumbline\Movie3_Scene1_PIV_00.raw'));
Combo_grid = IMgrid;
Combo_grid(:,5761:5800) = 5*double(imadjust(uint8(IMflatAfter2(:,5761:5800))));
Combo_grid(771:840,1:5000) = 4*double(imadjust(uint8(IMpl(771:840,1:5000)))); Combo_grid(5315:5360,1:4300) = 4*double(imadjust(uint8(IMpl(5315:5360,1:4300))));

% Compute Combo_grid undistorted
Combo_grid_lens = CorrectPIVLensDistortion(Combo_grid);
figure;imagesc(Combo_grid_lens);colormap gray

% Compute Combo_grid with camera angle correction to check straight lines,
% vertical plumblines and horizontal flat surface
Combo_grid_CamAngle_Corrected = Correct_Angle_PIV(Combo_grid_lens);
figure;imagesc(Combo_grid_CamAngle_Corrected)

%% IMgrid used from PIVSURF to match PIV resolution
IMgrid_lens = CorrectPIVLensDistortion(IMgrid);
IMgrid_CamAngle_Corrected = Correct_Angle_PIV(IMgrid_lens);

% Find checkerboard corners
[imagePoints, boardSize] = detectCheckerboardPoints(uint16(IMgrid_CamAngle_Corrected));
PIV_iP = imagePoints(32:end-31*12,:);
figure;imagesc(IMgrid_CamAngle_Corrected);hold on;plot(PIV_iP(:,1),PIV_iP(:,2),'ro')

% Points used: 1) Q-9 (O left-under) ; 2) B-9 (B left-under) ; 3) B-22 (B left-up) ; 4) Q-22 (O left-up)
Xp = [imagePoints(cursor_info(1).DataIndex,:) ; imagePoints(cursor_info(2).DataIndex,:);...
      imagePoints(cursor_info(3).DataIndex,:) ; imagePoints(cursor_info(4).DataIndex,:)];
  
% Save matching points
save('\\spray1-10g\d\Shoaling_waves\M-files\Step 3 - Rectify & Surface Detection\calibration\PIV_PIVsurf_matching_points', '-append', 'Xp', 'PIVgrid_CamAngle_Corrected', 'PIV_iP')