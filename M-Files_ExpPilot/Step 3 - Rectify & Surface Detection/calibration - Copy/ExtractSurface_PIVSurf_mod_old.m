function [XPos,YPos,PIVSurf_PIVMatch,XPIV_PIVSurf,PIVSurf_Surface,PIV_match_Surface,BadFramePIVSurf] = ExtractSurface_PIVSurf_mod(PIVSURF_lens,FusedPIV1)

%%%%%%%% This function extracts PIVSurf surface and matches PIV coordinates
%
%%% Input: 
%       LFV_lens = LFV image after lens distortion correction
%       XPos = XPos obtained from PIVSurf transformation to match PIV
%              coordinates
%       YPos = YPos obtained from PIVSurf transformation to match PIV
%              coordinates
%       XPIV_PIVSurf_Surface = X-coordinate in PIV resolution for PIVSurf
%                              surface
%       PIVSurf_Surface = PIVSurf surface position in PIV resolution
%       PIVFused_Surface = PIVSurf surface exactly matching PIV length
% 
%%% Output: 
%       LFV_lens = LFV image after lens distortion correction
%       XPos = XPos obtained from PIVSurf transformation to match PIV
%              coordinates
%       YPos = YPos obtained from PIVSurf transformation to match PIV
%              coordinates
% 

%% Fabrice's old PIVSurf (for URI2918 experiments)
% 
%%%%% Old Fabrice's function
% function [BadFramePIVSurf,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface,PIVSurf_PIVMatch] = ExtractSurface_PIVSurf(NORM_PIV_SURF,PIVSurf_CRR,FusedPIV1)
% 
% % % % Y1 = abs(round(PIVSurf_CRR.Ypos(1))); Y2 = Y1 + size(FusedPIV1,1) - 1;
% % % % X1 = abs(round(PIVSurf_CRR.Xpos(1))); X2 = X1 + size(FusedPIV1,2) - 1;
% % % % %Coordinate matching PIV image
% % % % 
% % % % h = fspecial('gaussian',16,10);
% % % % PIVSurf_CRR.img(:,4645:end)=PIVSurf_CRR.img(:,4645:end)-11; %offset correction for second tap (right side of image)
% % % % PIVSurf_PIVMatch=PIVSurf_CRR.img(Y1:Y2,X1:X2);
% % % % NS=NORM_PIV_SURF(:,219:end-256);
% % % % f=NS(1,:); NS=[repmat(f,1500,1); NS];
% % % % PIVSurf_PIVMatch=PIVSurf_PIVMatch./NS;
% % % % PIVSurf_PIVMatch=imfilter(PIVSurf_PIVMatch,h,'replicate');
% % % % %PIVSurf image exactely matching PIV images
% % % % 
% % % % %PIVSurf_Strip = PIVSurf_CRR.img(Y1+1500:Y2, 1+361:end-361);  %Correct size for matlab 2014
% % % % %%%ATTENTION%%% use the following instead if using Matlab version 2015 and
% % % % %%%above! NORM_PIV_SURF was obtained with a version 2014 and is larger by
% % % % %%%two columns - There was a change in "IMROTATE" routine use line 34 of
% % % % %%%CorrectPIVSurfv_v2.
% % % % %%NN=361 for old Matlab NN=360 for new one
% % % % NN=360;
% % % % PIVSurf_Strip = PIVSurf_CRR.img(Y1+1500:Y2, 1+NN:end-NN);  %Correct size for matlabe 2015 and later
% % % % PIVSurf_Strip=PIVSurf_Strip./NORM_PIV_SURF;
% % % % PIVSurf_Strip=imfilter(PIVSurf_Strip,h,'replicate');
% % % % 
% % % % % PIVSurf image exactely matching height of PIV images, but longer
% % % % % 361 gets rid of NaN from the rotation and matching of PIVSurf with PIV images
% % % % %1500 is a vertical offset to get closer to the surface and reduce image
% % % % %size in which to detect the surface
% % % % PF_XPos_in_PIVSurf = X1-360:X2-360; %to re-place surface from strip into PIV coordinates
% % % % 
% % % % %Creating filter around surface
% % % % PIVSurf_Strip_T = PIVSurf_CRR.img(Y1+1550:Y2+50, 1+NN:end-NN);
% % % % PIVSurf_Strip_B = PIVSurf_CRR.img(Y1+1450:Y2-50, 1+NN:end-NN);
% % % % PIVSurf_Strip_T=PIVSurf_Strip_T./NORM_PIV_SURF;
% % % % PIVSurf_Strip_B=PIVSurf_Strip_B./NORM_PIV_SURF;
% % % % S_T=(imfilter(PIVSurf_Strip_T,fspecial('gaussian',64,16),'replicate'));
% % % % S_B=(imfilter(PIVSurf_Strip_B,fspecial('gaussian',64,16),'replicate'));
% % % % S=((S_T-S_B));
% % % % 
% % % % PIVSurf_Surface_Tem = FindSurface(PIVSurf_Strip,5,S); %5
% % % % %finds surface in the strip
% % % % 
% % % % BadFramePIVSurf=PIVSurf_Surface_Tem.badFrameBool;
% % % % 
% % % % PIVSurf_Surface_Raw = PIVSurf_Surface_Tem.surface; 
% % % % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% % % % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% % % % PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
% % % %         
% % % %         try
% % % %             PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);       
% % % %         catch
% % % %             warning('on');
% % % %             warning(['Problem using FiltSpray function. Ignoring the ' ...
% % % %                      'FiltSpray function for PIVSurf image']);
% % % %             PIVSurf_Surface_Int = PIVSurf_Surface_Raw;
% % % %             BadFramePIVSurf=1;
% % % %         end
% % % %         
% % % %        PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
% % % %         % Interpolates NaN in the surface
% % % %                 
% % % %         [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 100000); %100000
% % % %         %Smoothing of the surface
% % % %         if (length(SP.coefs)>length(PIVSurf_Surface_Int))
% % % %             PIVSurf_Surface = SP.coefs(2:end-1);
% % % %         else
% % % %             PIVSurf_Surface=PIVSurf_Surface_Int;
% % % %         end
% % % %         PIVSurf_Surface=PIVSurf_Surface+1500; %Correction including vertical offset 
% % % %         PIVFused_Surface = PIVSurf_Surface(PF_XPos_in_PIVSurf); 
% % % %         % Portion of surface that corresponds to PIV 
% % % %         XPIV_PIVSurf_Surface=[-(X1-361-1):1:length(PIVSurf_Surface)-(X1-361-1)-1];
% % % %        
% % % %         
% % % % if(max(PIVFused_Surface) > size(FusedPIV1,1) || min(PIVFused_Surface) < 1)
% % % %     BadFramePIVSurf=1;
% % % % end

%% PIVSurf
% Polynomial transformation to correct for residual errors from camera
% distortion
load PIV_match.mat U1 X1 Xt
T1 = fitgeotrans(U1,X1,'polynomial',4);
PIVSURF_CamAngle = PIVSURF_lens;
PIVSURF_CamAngle(isnan(PIVSURF_CamAngle)) = 0;
PIVSURF_CamAngle = imwarp(PIVSURF_CamAngle,T1,'linear');
PIVSURF_CamAngle = imtranslate(PIVSURF_CamAngle, [(X1(1,1)-Xt(1,1)) (X1(1,2)-Xt(1,2))] );
PIVSURF_CamAngle = fliplr(PIVSURF_CamAngle);

% Find transformation to match pixel's resolution on PIV image
P = load('PIV_PIVsurf_matching_points.mat','Xp','Up','PIV_iP','PIVSURF_iP','PIVSURF_iP2');
Up2 = [-X1(:,1)+size(PIVSURF_CamAngle,2)+1,X1(:,2)];
Up2([1:26,end-8*26+1:end],:) = [];
Up2 = flipud(reshape(permute(reshape(Up2,[26,31,2]),[2,1,3]),26*31,2));
Up = [ Up2(end-30,:) ; Up2(end,:) ; Up2(31,:) ; Up2(1,:) ]; % The coordinate
% of the points in the PIVSurf image (the grid calibration image).
Xp = P.Xp; % The coordinates of the points on PIV image. The physical location  of these points are
% the same as PIVSurf images (Up). The calibration grid was used to find out
% the exact same locations for these two images.
T2_bis = maketform('projective',Up,Xp);

% PIVSurf image matching PIV
[Resized_PIVSurf_bis,XPos,YPos] =  imtransform(PIVSURF_CamAngle,T2_bis,'XYScale',1);
PIVSurf_PIVMatch = Resized_PIVSurf_bis(1:size(FusedPIV1,1)-YPos(1),1-XPos(1):size(FusedPIV1,2)-XPos(1));
PIVSurf_PIVMatch = [zeros(round(YPos(1)),size(PIVSurf_PIVMatch,2));PIVSurf_PIVMatch];

%% Find surface on resized PIVSurf

% Step 1: find surface on PIVSurf image with original size
X = 1:size(PIVSURF_CamAngle,2);
Xgood = 301:length(X); % good points where image is still bright enough to detect surface
X = Xgood;
[imSurf_PIVSURF] = FindSurface(PIVSURF_CamAngle,5,1);
PIVSurf_Surface_Raw = despike_fab(imSurf_PIVSURF.surface(X));
PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1000); %100000
% PIVSurf_Surface = SP.coefs; %%% use for flat surface
PIVSurf_Surface_2 = SP.coefs(2:end-1);
% PIVSurf_Surface = polyval(polyfit([X(1),X(end)],PIVSurf_Surface,1),X);  %%% use for flat surface
PIVSurf_Surface_resized2 = tformfwd([X',PIVSurf_Surface_2'],T2_bis);
%         PIVSurf_Surface_resized = interp1(PIVSurf_Surface_resized2(:,1),PIVSurf_Surface_resized2(:,2),1:length(PIV_Surface));
X_PIVSURF2 = PIVSurf_Surface_resized2(:,1);
PIVSurf_Surface_resized = PIVSurf_Surface_resized2(:,2);
PIVSurf_SurfMatch_PIV = interp1(PIVSurf_Surface_resized2(:,1),PIVSurf_Surface_resized2(:,2),1:size(FusedPIV1,2));

% Further rotation to correct flat surface
xrr = PIVSurf_Surface_resized';
XX = [X_PIVSURF2';xrr];
Alpha = rad2deg(9.15/5924);
DX = 25.3894;%34.5595;
Yref = 5626.3;
MM = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
Xr = MM*XX;
%         xr{i+1} = interp1(Xr(1,:),Xr(2,:),1:5924,'pchip','extrap');
XPIV_PIVSurf = round(X_PIVSURF2(1)):round(X_PIVSURF2(end));
PIVSurf_Surface = interp1(Xr(1,:),Xr(2,:),XPIV_PIVSurf,'pchip','extrap')-DX;

% PIVSurf surface matching PIV coordinates
IndX = [find(XPIV_PIVSurf==1),find(XPIV_PIVSurf==size(FusedPIV1,2))];
PIV_match_Surface = PIVSurf_Surface(IndX(1):IndX(end));

%% Check bad frames
BadFramePIVSurf = [];
if (max(PIV_match_Surface) > size(FusedPIV1,1) || min(PIV_match_Surface) < 1)
    BadFramePIVSurf = 1;
end