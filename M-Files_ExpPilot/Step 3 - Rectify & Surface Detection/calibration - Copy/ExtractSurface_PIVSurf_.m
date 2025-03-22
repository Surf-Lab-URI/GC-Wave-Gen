function [BadFramePIVSurf,PIVSurf_PIV,X,PIVSurf_Surf,X_PIVSURF,PIVSURF_PIVMatch,IndX] = ExtractSurface_PIVSurf_(PIVSURF_CamAngle,PIV1)

%% PIVSURF surface
% (1) Extract surface
% (2) Match PIV coordinates

%% Correct camera angle
load PIV_match.mat X1 

%% Resize PIVSURF on PIV image
% Resize PIV image and match PIV coordinates
[PIVSurf_Corrected,T2_bis] = Resize_image_PIVSURF(PIVSURF_CamAngle);

% Crop PIVSURF to exactly match PIV image
PIVSurf_PIV = PIVSurf_Corrected.img(1:size(PIV1,1)-PIVSurf_Corrected.YPos(1),1-PIVSurf_Corrected.XPos(1):size(PIV1,2)-PIVSurf_Corrected.XPos(1));
PIVSurf_PIV = [zeros(round(PIVSurf_Corrected.YPos(1)),size(PIVSurf_PIV,2));PIVSurf_PIV];

%% Find surface on original PIVSURF
X = 1:size(PIVSURF_CamAngle,2);
Xgood = 301:length(X); % interval where laser light is still strong enough to illuminate surface
X = Xgood;

% Find surface on PIVSURF of original size (after correcting lens distortion and angle)
[imSurf_PIVSURF] = FindSurface(PIVSURF_CamAngle,5,1);
PIVSurf_Surf_Raw = despike_fab(imSurf_PIVSURF.surface(X));
PIVSurf_Surf_Raw = despike_fab(PIVSurf_Surf_Raw);
PIVSurf_Surf_Raw = despike_fab(PIVSurf_Surf_Raw);
PIVSurf_Surf_Int = filt_spray(PIVSurf_Surf_Raw);
PIVSurf_Surf_Int = smoothn(PIVSurf_Surf_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surf_Int), PIVSurf_Surf_Int, 1000); %100000

% Surface on PIVSURF size 
PIVSurf_Surf = SP.coefs(2:end-1);

%%% Use ONLY for flat surface
%PIVSurf_Surf = SP.coefs;
%PIVSurf_Surf = polyval(polyfit([X(1),X(end)],PIVSurf_Surf,1),X);

%% Transform surface to match resized PIVSURF
% Resize PIVSURF surface on PIV
PIVSurf_Surface_resized2 = tformfwd([X',PIVSurf_Surf'],T2_bis);
X_PIVSURF2 = PIVSurf_Surface_resized2(:,1);
PIVSurf_Surface_resized = PIVSurf_Surface_resized2(:,2);

% PIVSURF surface matching PIV coordinates
% % % PIVSurf_SurfMatch_PIV = interp1(PIVSurf_Surface_resized2(:,1),PIVSurf_Surface_resized2(:,2),1:length(PIV_Surface));

%%% Check on PIVSURF image
% figure
% imagesc()
% hold on
% plot()

%% Further rotation to correct flat surface (matched for 8cm w.d.)
xrr = PIVSurf_Surface_resized';
XX = [X_PIVSURF2';xrr];
Alpha = rad2deg(9.15/5924);
DX = 25.3894;
Yref = 5626.3;
MM = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
Xr = MM*XX;
X_PIVSURF = round(X_PIVSURF2(1)):round(X_PIVSURF2(end));
PIVSURF_PIVMatch = interp1(Xr(1,:),Xr(2,:),X_PIVSURF,'pchip','extrap')-DX;
IndX = [find(X_PIVSURF==1),find(X_PIVSURF==size(PIV1,2))]; %%% CHECK IF THE 1ST COLUMN IS 0 OR 1!!!

%% Check bad frames
BadFramePIVSurf = [];
if(max(PIVSURF_PIVMatch) > size(PIV1,1) || min(PIVSURF_PIVMatch) < 1)
    BadFramePIVSurf = 1;
end