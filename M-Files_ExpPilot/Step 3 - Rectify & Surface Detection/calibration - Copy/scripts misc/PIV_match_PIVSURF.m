%%% Finding calibration points for dewarping of PIV
clear 
close all 
clc

% Load image and detect checkerboard points
% % % [IM2] =
% fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie2_Scene4\RAW\PIV\Movie2_Scene4_PIV_0.raw')); %old
IM2 = fliplr(load_Image_IOCoreView_12MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\PIV_match\05232022_plumbline_Movie1_Scene3_PIVSURF_0.raw'));
[IM] = CorrectPIVLensDistortion(IM2);
IM(isnan(IM)) = 0;
[imagePoints, boardSize] = detectCheckerboardPoints(uint16(IM));
figure;imagesc(IM);colormap gray;hold on;plot(imagePoints(:,1),imagePoints(:,2),'ro')

% Reshape points and find mean distance
% % % IMpoints = reshape(imagePoints,35,39,2); old
Col = boardSize(1)-2;
Row = boardSize(2)-1;
imagePoints(Col+1:Col+1:end,:) = []; 
IMpoints = reshape(imagePoints,Col,Row,2);
diffIMpoints = diff(IMpoints);

% Show trend in the vertical and horizontal
figure;hold on
for i = 1:Col %35
    plot(diff(IMpoints(i,:,1)))
end
figure;hold on
for i = 1:Row %39
    plot(diff(IMpoints(:,i,2)))
end

% Compute mean distance in the vertical and in the horizontal
for i =1:Col %35
    dxM(i) = nanmean(diff(IMpoints(i,:,1)));
end
for i =1:Row %39
    dyM(i) = nanmean(diff(IMpoints(:,i,2)));
end

% Show mean distance in the vertical and in the horizontal 
figure;plot(dxM)
figure;plot(dyM)

% Find mean reprojection distance for both vertical and horizontal
Dpix = mean([mean(dxM),mean(dyM)]);
Dstd = [mean(dxM),mean(dyM);std(dxM),std(dyM)];

%%% Construct points for fitgeotrans
% Detected points
U1 = imagePoints;

% Reprojected points
X1 = nan(size(U1,1),size(U1,2));
for i = Row:-1:1 %39
    X1(1+(i-1)*Col,1) = U1(end,1)-Dpix*(i-1);
    X1(1+(i-1)*Col:i*Col,1) = X1(1+(i-1)*Col,1);
end
for i = Col:-1:1 %35
    X1(i,2) = U1(end,2)-Dpix*(i-1);
    X1(i:Col:end,2) = X1(i,2);
end

figure;imagesc(IM);colormap gray;hold on;plot(U1(:,1),U1(:,2),'ro',X1(:,1),X1(:,2),'bx')

% Save results
save('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIVSURF\PIV_match\PIV_match','U1','X1','imagePoints','IM','IMpoints','dxM','dyM','Dpix');