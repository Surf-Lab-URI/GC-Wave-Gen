clear 
close all
clc

% Find flat surface on the right and on the left
IMgrid = double(imread('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\PIVSURF_match\Movie2_Scene4_PIV_2.png'));
IMflat = fliplr(load_Image_IOCoreView_48MP('\\spray1-10g\d\Shoaling_waves\data\Calibration\PIV\flat\Movie1_Scene1_PIV_000.raw'));
IM_lens = CorrectPIVLensDistortion(IMflat);
[IM_CamAngle_Corrected] = Correct_Angle_PIV(IM_lens);
img = double(imadjust(uint8(IM_CamAngle_Corrected)));
figure;imagesc(img)
x = nanmean(img(5901:end,:));
y = nanmean(img(51:150,:));
[r,lag] = xcorr(x,y);
[~,Imax] = max(r);
[lag(Imax)]
figure;plot(lag,r,lag(Imax),r(Imax),'ro')
img(img==0)=nan;
[imSurf] = FindSurface(img',5,1);
imSurf.surface(1:150) = nan;
figure;imagesc(img);hold on;plot(movavg(interp1(1:length(imSurf.surface),imSurf.surface,1:length(imSurf.surface),'pchip','extrap')','linear',2120),1:length(imSurf.surface),'r')

%%% Find mean gaussian peak for all the plumbline images in the final image
%%% for further projective transformation
for i = 2:5
IM = fliplr(load_Image_IOCoreView_48MP(['\\spray1-10g\d\Shoaling_waves\data\Calibration\05182022_flat\Movie1_Scene' num2str(i) '\RAW\PIV\Movie1_Scene' num2str(i) '_PIV_00.raw']));
IM_lens = CorrectPIVLensDistortion(IM);
[IM_CamAngle_Corrected] = Correct_Angle_PIV(IM_lens);
x_T = nanmean(IM_CamAngle_Corrected(:,101:200),2);
x_B = nanmean(IM_CamAngle_Corrected(:,3941:4040),2);
T = findpeaks(x_T); T = T(T>400); PT_R(i) = T(1); PT_L(i) = T(2);
B = findpeaks(x_B); B = B(B>400); PB_R(i) = B(1); PB_L(i) = B(2);
keyboard
end

%%% Find mean gaussian peak for all the plumbline images before the angle
%%% correction
for i = 2:5
IM = fliplr(load_Image_IOCoreView_48MP(['\\spray1-10g\d\Shoaling_waves\data\Calibration\05182022_flat\Movie1_Scene' num2str(i) '\RAW\PIV\Movie1_Scene' num2str(i) '_PIV_00.raw']));
Combo_flat = IM;
Combo_flat(:,5601:5900) = 4*double(imadjust(uint8(IMflat(:,5601:5900))));
Combo_flat_lens = CorrectPIVLensDistortion(Combo_flat);
x_T = nanmean(Combo_flat_lens(:,101:200),2);
x_B = nanmean(Combo_flat_lens(:,3941:4040),2);
figure;plot(nanmean(Combo_flat_lens(:,101:200),2),'r');hold on;plot(nanmean(Combo_flat_lens(:,3941:4040),2),'b')
keyboard
end