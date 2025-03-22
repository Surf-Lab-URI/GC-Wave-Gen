%% Import RAW calibration images and write FMT (png, jpg, tif) images
% This script was used to write images that can be uploaded by
% cameraCalibrator App (or by Caltech toolbox)

clear
close all
clc

%%%% Lens Calibration bigger grid LFV
% % Import folder
% CalDir = dir('Y:\Shoaling_waves\data\Calibration\LFV\chessboard\*.raw');
% IM = cell(length(CalDir),1);
% for i = 1:length(CalDir)
%     IM{i} = CalDir(i).name;
% end
% 
% for IndIm = 1:length(CalDir)
%     eval(['[IM' num2str(IndIm) '] = load_Image_IOCoreView([CalDir(IndIm).folder ''\'' IM{IndIm}]);']);
%     %figure;eval(['imagesc(IM' num2str(IndIm) ')']);axis on;colormap gray
%     ImName = [CalDir(IndIm).folder '\' IM{IndIm}(1:end-4) '.png'];
%     eval(['imwrite(uint16(IM' num2str(IndIm) '*4),ImName, ''png'')'])
% end
%%%%

%Images directory
% CalDir = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie1_Scene3\RAW\LFV\*.raw');
% CalDir = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie1_Scene2\RAW\PIVSURF\*.raw');
% CalDir = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie1_Scene1\RAW\PIV\*.raw');
% CalDir = dir('\\spray1\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie3_Scene1\RAW\PIV\*.raw');
% CalDir = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05162022_calibration\Movie2_Scene2\RAW\PIV\*.raw');
CalDir1 = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05232022_calibration\Movie1_Scene2\RAW\PIV\*.raw');
CalDir2 = dir('\\spray1-10g\d\Shoaling_waves\data\Calibration\05232022_calibration\Movie1_Scene4\RAW\PIV\*.raw');

% Folder where to save images
% SaveFold = '\\spray1-10g\d\Shoaling_waves\data\Calibration\LFV\Checker board smaller grid\';
% SaveFold = '\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\checkerboard\';
% SaveFold = '\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\checkerboard\smaller\';
% SaveFold = '\\spray1\d\Shoaling_waves\data\Calibration\PIV\PIVSURF_match\';
% SaveFold = '\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\checkerboard\05162022\';
SaveFold = '\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\checkerboard\05232022\';

% Number of images to extract
% Int = [299,335,350,370,395,424,448,453,488,508,526,590,618,678,688,707,726,740,755,761,780,806,821,842,853,872,885,898,917,935,975,1000,1015,1040,1055,1071,1081,1093,1104,1116,1124,1140,1150,1160,1270]; % calibration LFV smaller grid 
% Int = [996,1053,1100,1141,1175,1178,1195,1250,1290,1330,1342,1361,1380,1383,1400,1420,1425,1434,1468,1475,1493,1515,1543,1570,1628,1678,1700,1709,1747,1785,1805,1816,1843,1865,1876,1882,1928,1988,2007,2010,2042,2075,2111,2188,2215,2225,2273,2321,2343,2381,2404,2430,2447,2477,2521,2553,2584,2612,2659,2697,2726,2764]; % calibration PIVSURF
% Int = [603,615,625,629,638,639,673,777,833,864,870,934,950,996,1009,1015,1150,1175,1197,1213,1225,1233,1247,1260,1376,1382,1469,1475,1506,1511,1596,1608,1613,1629,1635,1662,1861,1868,1875,1885,1904,1927,1933]; % calibration PIVSURF
% Int = [1]; % PIV-PIVSURF match grid
% Int1 = [419 425 447 449 463 482 510 578 582 588 609 647 657 685 706 717 868 935 961 994 1023 1027 1048 1055 1137 1142 1152 1156 1264 1268 1274 ]; % PIV lens dist 05162022_calibration
Int1 = [332,335,367,378,400,429,453,494,538,562,573,612,635,706,718,795,803,832,875,896,917,933,967,1027,1037,1126,1159,1227,1248,1262]; % PIV lens dist 1st-part 05232022_calibration
Int2 = [303,314,326,360,370,391,401,412,452,484,522,531,592,607,627,668,695,731,776,801,871,875,903,950,982,1007,1042]; % PIV lens dist 2nd-part 05232022_calibration

% Camera Megapixels 
MP = 48;

% Format of the image to write
FMT = 'tif';
% FMT = 'jpg';

%% Write Images from Raw
calibration_RAWtoFMT(CalDir1,SaveFold,Int1,MP,FMT);
calibration_RAWtoFMT(CalDir2,SaveFold,Int2,MP,FMT);