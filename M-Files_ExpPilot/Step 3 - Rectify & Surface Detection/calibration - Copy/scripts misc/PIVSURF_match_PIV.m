%%% Finding calibration points for dewarping of PIV
clear 
close all 
clc

% Load image and detect checkerboard points
% % % [IM2] =
% fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\05112022_calibration\Movie2_Scene4\RAW\PIV\Movie2_Scene4_PIV_0.raw')); %old
[IM2] = fliplr(load_Image_IOCoreView_48MP('\\spray1\d\Shoaling_waves\data\Calibration\PIV\PIVSURF_match\Movie1_Scene1_PIV_00.raw'));
[IM] = CorrectPIVLensDistortion(IM2);
IM(isnan(IM)) = 0;
[imagePoints, boardSize] = detectCheckerboardPoints(uint16(IM));
figure;imagesc(IM);colormap gray;hold on;plot(imagePoints(:,1),imagePoints(:,2),'ro')

% Reshape points and find mean distance
% % % IMpoints = reshape(imagePoints,35,39,2); old
Col = boardSize(1)-1;
Row = boardSize(2)-2;
imagePoints = imagePoints(1:Col*Row,:);
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

%% Construct points for fitgeotrans using fixed square dimensions
% Detected points
U1 = flipud(imagePoints);

% Reprojected points
X1 = nan(size(U1,1),size(U1,2));
for i = 1:Row %39
    X1(1+(i-1)*Col,1) = U1(1,1)-Dpix*(i-1);
    X1(1+(i-1)*Col:i*Col,1) = X1(1+(i-1)*Col,1);
end
for i = 1:Col %35
    X1(i,2) = U1(1,2)-Dpix*(i-1);
    X1(i:Col:end,2) = X1(i,2);
end

figure;imagesc(IM);colormap gray;hold on;plot(U1(:,1),U1(:,2),'ro',X1(:,1),X1(:,2),'bx')

%% Using interpolation to correct distortion (avoiding guessing grid square size and then reducing the errors)
cc1 = [];
cc2 = [];
% Correct residual distortion
for i = 1:31
    cc1{i} = polyfit(IMpoints(i,:,1),IMpoints(i,:,2),1);
end
for i = 1:39
    cc2{i} = polyfit(IMpoints(:,i,1),IMpoints(:,i,2),1);
end
figure;imagesc(IM');hold on;plot(imagePoints(:,2),imagePoints(:,1),'ro')
for i = 1:39
    plot(polyval(cc2{i},1:7920),1:7920,'c','LineWidth',2)
end
for i = 1:31
    plot(polyval(cc1{i},1:7920),1:7920,'w','LineWidth',2)
end

for i = 1:31
    x1 = [1;7920];
    y1 = [polyval(cc1{i},1);polyval(cc1{i},7920)];
    for ii = 1:39
        x2 = [1;7920];
        y2 = [polyval(cc2{ii},1);polyval(cc2{ii},7920)];
        [xi(i,ii),yi(i,ii)] = polyxpoly(x1,y1,x2,y2);
    end
end
Reprojected_points = [reshape(yi,31*39,1),reshape(xi,31*39,1)];

%% Save results
save('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\PIVSURF_match\PIVSURF_match','U1','X1','imagePoints','IM','IMpoints','dxM','dyM','Dpix','Reprojected_points');