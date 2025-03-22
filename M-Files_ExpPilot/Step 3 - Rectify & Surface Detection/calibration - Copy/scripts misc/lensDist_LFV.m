% Lens distortion correction for LFV

clear
close all 
clc

% Calibration images folder
CalDir = dir('Y:\Shoaling_waves\data\Calibration\LFV\chessboard\*.png');
CalIm = cell(length(CalDir),1);
for i =1:length(CalDir)
    CalIm{i} = [CalDir(i).folder '\' CalDir(i).name];
end

%% Checkerboard points detection

% Detection of the 
[imagePoints2,boardSize] = detectCheckerboardPoints(CalIm);

% Flipping points taken in the wrong direction (hard-coded from Camera
% Calibrator)
imagePoints = imagePoints2;
indFlip = [1,2,3,4,5,10,11];
for i = indFlip
    imagePoints(:,:,i) = flipud(imagePoints2(:,:,i));
end

% Indexes of good points (hard-coded checking from Camera Calibrator)
IndGood = [2,3,5,6,7,16,17,18,21,30,31,34,35];

% Checking the consistency of the detected points
for i = 1:size(imagePoints,3)
    figure
    Im = imread(CalIm{IndGood(i)});
    imagesc(Im)
    colormap gray
    hold on
    plot(imagePoints(:,1,i),imagePoints(:,2,i),'ro')
end

% World coordinates for the detected points
squareSize = 2*unitsratio('mm','in');
worldPoints = generateCheckerboardPoints(boardSize,squareSize);

%% Estimate parameters

imageSize = [size(Im, 1),size(Im, 2)];
[params,~,errors] = estimateCameraParameters(imagePoints,worldPoints,'ImageSize',imageSize,'NumRadialDistortionCoefficients',3);

% Visualize the calibration accuracy
figure;
showReprojectionErrors(params);

%% Save parameters
save([CalDir(1).folder,'\params.mat'],'params','errors')