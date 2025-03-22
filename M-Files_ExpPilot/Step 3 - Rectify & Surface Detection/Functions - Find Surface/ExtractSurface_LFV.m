function [BadFrameLFVSurf, XPIV_LFV,LFV_Surface_PIVMatched,LFV_Surface_Combo_Surface] = ExtractSurface_LFV(LFV_CR1,XPos,YPos,PIV1,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIV_match_Surface,expName)

%%%%%%%% This function extracts LFV surface and matches PIV coordinates
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

%% Fabrice's old LFV (for URI2918 experiments)

% % Old function used by Fabrice
% function [BadFrameLFVSurf, XPIV_LFV_Surface_PIVMatched,LFV_Surface_PIVMatched,LFV_Surface_Combo_Surface] = ExtractSurface_LFV(NORM_LFV,LFV_CR1,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface)
% 

% % % % % % LFV_CR2=LFV_CR1./NORM_LFV;
% % % % % % h = fspecial('gaussian',16,10);
% % % % % % LFV_CR3=imfilter(LFV_CR2,h,'replicate');
% % % % % % 
% % % % % % LFV_CR4=LFV_CR3(200:1800,157:1680); %to get rid of NAN from lens distortion correction
% % % % % % 
% % % % % % S_T=(imfilter(LFV_CR2(250:1850,157:1680),fspecial('gaussian',64,16),'replicate'));
% % % % % % S_B=(imfilter(LFV_CR2(150:1750,157:1680),fspecial('gaussian',64,16),'replicate'));
% % % % % % S=((S_T-S_B));
% % % % % % 
% % % % % % LFV_Surface_Tem = FindSurface(LFV_CR4, 1,S);
% % % % % % %finds surface in the strip
% % % % % % 
% % % % % % BadFrameLFVSurf=LFV_Surface_Tem.badFrameBool;
% % % % % % 
% % % % % % LFV_Surface_Raw = LFV_Surface_Tem.surface; 
% % % % % % LFV_Surface_Raw=despike_fab(LFV_Surface_Raw);
% % % % % % LFV_Surface_Raw=despike_fab(LFV_Surface_Raw);
% % % % % % LFV_Surface_Raw=despike_fab(LFV_Surface_Raw);
% % % % % %         
% % % % % %         try
% % % % % %             LFV_Surface_Int = filt_spray(LFV_Surface_Raw);       
% % % % % %         catch
% % % % % %             warning('on');
% % % % % %             warning(['Problem using FiltSpray function. Ignoring the ' ...
% % % % % %                      'FiltSpray function for PIVSurf image ']);
% % % % % %             LFV_Surface_Int = LFV_Surface_Raw;
% % % % % %             BadFrameLFVSurf=1;
% % % % % %         end
% % % % % %         
% % % % % %        LFV_Surface_Int=LFV_Surface_Int(1:length(LFV_Surface_Raw));
% % % % % %        
% % % % % %        LFV_Surface_Int = smoothn(LFV_Surface_Int, 'robust');
% % % % % %         % Interpolates NaN in the surface
% % % % % %                 
% % % % % %         [SP,~] = spaps(1:length(LFV_Surface_Int), LFV_Surface_Int, 400);
% % % % % %         %Smoothing of the surface
% % % % % %         if (length(SP.coefs)>length(LFV_Surface_Int))
% % % % % %             LFV_Surface = SP.coefs(2:end-1);
% % % % % %         else
% % % % % %             LFV_Surface=LFV_Surface_Int;
% % % % % %         end
% % % % % %         
% % % % % %         if(max(LFV_Surface) > size(S,1) || min(LFV_Surface) < 1)
% % % % % %        BadFrameLFVSurf=1;
% % % % % %         end
% % % % % % 
% % % % % %         
% % % % % %        load LFV_flat_correction.dat
% % % % % %        X_LFV=[1:size(LFV_CR1,2)];
% % % % % %        LFV_Surface=[[1:156]*nan LFV_Surface+199 [1681:size(LFV_CR1,2)]*nan];
% % % % % %        LFV_Surface=LFV_Surface-LFV_flat_correction;
% % % % % %         % lens Correction again (from line 7), 199 accounts for the 200:1800 strip
% % % % % %         % 156 and 1681 acount for the 157:1680 strip
% % % % % %          
% % % % % %         load LFV_x_inches.dat
% % % % % %         %Contains x in inches (starting from 0) for the 2048 points LFV surface
% % % % % %         %That's used because of the lens distortion still present (parallax)
% % % % % %         XPIV_LFV=(LFV_x_inches-LFV_x_inches(1038))*590.3+3993;%point 1038 in LVF surface is actually point 3993 in PIV
% % % % % %         %XPIV_LFV is 2048 points long but is the location of the LFVsurface in PIV (in PIV points)
% % % % % %         XPIV_LFV_Surface_PIVMatched=[round(min(XPIV_LFV)):round(max(XPIV_LFV))];
% % % % % %         LFV_Surface_PIVMatched = interp1(XPIV_LFV,LFV_Surface,XPIV_LFV_Surface_PIVMatched,'linear','extrap');
% % % % % %         LFV_Surface_PIVMatched=(LFV_Surface_PIVMatched-1081)*590.3/58+2438;
% % % % % %         %1081 is flat level in LFV images
% % % % % %         %2438 is flat level in PIV images
% % % % % %         %590.3/58 is ratio of resolutions; i.e. 590.3 pix=1inch in PIV
% % % % % %         %image and 58pix=1inch in LFV (vertical AND near image center)
% % % % % %      %%making combo surface
% % % % % % index_LL=find(XPIV_LFV_Surface_PIVMatched==XPIV_PIVSurf_Surface(1));
% % % % % % index_L=find(XPIV_LFV_Surface_PIVMatched==0);
% % % % % % index_R=find(XPIV_LFV_Surface_PIVMatched==length(PIVFused_Surface)+1);
% % % % % % index_RR=find(XPIV_LFV_Surface_PIVMatched==XPIV_PIVSurf_Surface(end));
% % % % % % 
% % % % % % a = LFV_Surface_PIVMatched(1:index_LL-1);
% % % % % % b_L=LFV_Surface_PIVMatched(index_LL-1:index_L);
% % % % % % b = LFV_Surface_PIVMatched(index_LL:index_RR);
% % % % % % b_R=LFV_Surface_PIVMatched(index_R:index_RR);
% % % % % % c = LFV_Surface_PIVMatched(index_RR+1:end);
% % % % % % PS_Mod = PIVSurf_Surface;
% % % % % % PS_Mod_L=PS_Mod(1:length(b_L));
% % % % % % PS_Mod_R=PS_Mod(length(PS_Mod)-length(b_R)+1:end);
% % % % % % 
% % % % % % a=a+(mean(PS_Mod_L-b_L)); b_L=b_L+(mean(PS_Mod_L-b_L));
% % % % % % c=c+(mean(PS_Mod_R-b_R)); b_R=b_R+(mean(PS_Mod_R-b_R));
% % % % % % 
% % % % % % TZ_L = length(b_L);
% % % % % % TZ_R = length(b_R);
% % % % % % PS_Mod(1:TZ_L) = 1/(TZ_L-1) * ((0:TZ_L-1) .* PIVSurf_Surface(1:TZ_L) + (TZ_L-1:-1:0) .* b_L(1:TZ_L));
% % % % % % PS_Mod(end-TZ_R+1:end) = 1/(TZ_R-1) * ((0:TZ_R-1) .* b_R(end-TZ_R+1:end) + (TZ_R-1:-1:0) .* PIVSurf_Surface(end-TZ_R+1:end));
% % % % % %         
% % % % % % LFV_Surface_Combo_Surface = [a PS_Mod c];
% % % % % % i=isnan(LFV_Surface_Combo_Surface);
% % % % % % LFV_Surface_Combo_Surface=despike_fab(LFV_Surface_Combo_Surface);
% % % % % % LFV_Surface_Combo_Surface(i)=NaN;
% % % % % % LFV_Surface_Combo_Surface(index_L+1:index_R-1)=PIVFused_Surface;
% % % % % % 
% % % % % % i=~isnan(LFV_Surface_Combo_Surface);
% % % % % % LFV_Surface_Combo_Surface=LFV_Surface_Combo_Surface(i);
% % % % % % XPIV_LFV_Surface_PIVMatched=XPIV_LFV_Surface_PIVMatched(i);
% % % % % % LFV_Surface_PIVMatched=LFV_Surface_PIVMatched(i);

%% LFV 

% Use undistorted image
LFV_CamAngle = LFV_CR1;
LFV_CamAngle(isnan(LFV_CamAngle)) = 0;

% Correct for straight plumb lines and flat surface
UU1 = [ 184.5  1457 ; 3788.5 1497 ; 3991.5 100 ; 1.75 100];
XX1 = [ (184.5+1.75)/2 1477 ; (3788.5+3991.5)/2 1477 ; (3788.5+3991.5)/2 100 ; (184.5+1.75)/2 100 ];
TT1 = fitgeotrans(UU1,XX1,'projective');
PIV2_CamAngle_Corrected = imwarp(LFV_CamAngle,TT1,'linear');

% Find transformation to match PIVSurf size (in pixel resolution)
Up22 = [ 965 1026 ; 2947 876.5 ; 3459 1419 ; 968 1417 ];
Xp22 = [ 1196 1887 ; 7951 1351 ; 9759 3182 ; 1190 3187 ];
T2_bis2 = maketform('projective',Up22,Xp22);

% LFV image matching PIVSurf
% % [Resized_LFV_bis,XPos2,YPos2] =  imtransform(PIV2_CamAngle_Corrected,T2_bis2,'XYScale',1);
% % LFV_PIVSurf = Resized_LFV_bis(1-YPos2(1):7242-YPos2(1),1-XPos2(1):9937-XPos2(1));
       
%% Find surface on resized LFV

% Step 1: find surface on undistorted LFV
%%% S filter computed on PIVSurf_PIVmatch image
F1 = 750; %5000;
F2 = 3000;
if str2double(expName)>24 && str2double(expName) < 37
    F2 = 3300;
end
if str2double(expName)>48 && str2double(expName) < 72
    F2 = 2300;
end
Lint = 100;
S_T = PIV2_CamAngle_Corrected(F1+Lint+1:F2,:); S_T = imfilter(S_T,fspecial('gaussian',64,16),'replicate');
S_B = PIV2_CamAngle_Corrected(F1+1:F2-Lint,:); S_B = imfilter(S_B,fspecial('gaussian',64,16),'replicate');
SS = S_T-S_B;
S_pos2 = SS;S_pos2(S_pos2<0)=0;
S_neg2 = SS;S_neg2(S_neg2>0)=0;
S = [zeros(F1+Lint,size(S_neg2,2));S_pos2;zeros(size(PIV2_CamAngle_Corrected,1)-F2,size(S_neg2,2))];
% S(1951:end,:) = 0;

%Y = (1:size(PIV2_CamAngle_Corrected,2));
Ygood = 801:3900; % interval where image is bright enough to retrieve the surface
if str2double(expName)>48 && str2double(expName)<55
    Ygood = 601:3600; % interval where image is bright enough to retrieve the surface
end
if str2double(expName)>56 && str2double(expName)<60
    Ygood = 301:4900; % interval where image is bright enough to retrieve the surface
end
Y = Ygood;
[imSurf_LFV] = FindSurface_LFV(S,5,1,expName);
LFV_Surface_Raw = despike_fab(imSurf_LFV.surface(Y));
LFV_Surface_Raw = despike_fab(LFV_Surface_Raw);
LFV_Surface_Raw = despike_fab(LFV_Surface_Raw);
LFV_Surface_Int = filt_spray(LFV_Surface_Raw);
LFV_Surface_Int = smoothn(LFV_Surface_Int, 'robust');
[SP,~] = spaps(1:length(LFV_Surface_Int), LFV_Surface_Int, 25000); % smoothing function: before it was 1000
if str2double(expName) > 36
    [SP,~] = spaps(1:length(LFV_Surface_Int), LFV_Surface_Int, 500000); % smoothing function: before it was 1000
end
LFV_Surface = SP.coefs(2:end-1);

% Step 2: resize to match PIV resolution
LFV_Surface_resized_PIVSURF2 = tformfwd([Y',LFV_Surface'],T2_bis2); % transformation to match PIVSurf resolution
X_LFV2 = LFV_Surface_resized_PIVSURF2(:,1)+XPos(1);
LFV_Surface_resized2 = LFV_Surface_resized_PIVSURF2(:,2)+YPos(1);
% % % LFV_Surface_resized_PIVSURF = interp1(LFV_Surface_resized_PIVSURF2(:,1),LFV_Surface_resized_PIVSURF2(:,2),1:9937);
% % % LFV_SurfMatch_PIV = LFV_Surface_resized_PIVSURF(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1))+YPos(1);
XPIV_LFV = -1339:7852; % hard-coded to have all the Combo surfaces of the same length for all the images pairs; it derives from "round(X_LFV2(1)):round(X_LFV2(end)); "
if str2double(expName)>48 && str2double(expName)<61
    XPIV_LFV = -2000:7999; % hard-coded to have all the Combo surfaces of the same length for all the images pairs; it derives from "round(X_LFV2(1)):round(X_LFV2(end)); "
end
LFV_Surface_resized = interp1(X_LFV2,LFV_Surface_resized2,XPIV_LFV);
IndX2 = [find(XPIV_LFV==1),find(XPIV_LFV==size(PIV1,2))];

%%% Step 3: further rotation to correct flat surface
xrr = LFV_Surface_resized;
XX = [(XPIV_LFV);xrr];
Alpha = rad2deg(9.15/5924);
DX2 = 30.3151; % vertical translation to match real position (based on empirical observation of random images)
MM = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
Xr = MM*XX;
Xr(:,isnan(Xr(1,:))) = [];
LFV_Surface_PIVMatched = interp1(Xr(1,:),Xr(2,:),XPIV_LFV,'pchip','extrap'); % final LFV surface with PIV image resolution
DDX = 57.4641;
DDDX = 0;
if str2double(expName)>48
    DDX = 0;
    DDDX = 470;
end
LFV_Surface_PIVMatched = LFV_Surface_PIVMatched-DDX;
XPIV_LFV = XPIV_LFV-DDDX;

%% Making combo surface with PIVSurf
index_LL=find(XPIV_LFV==XPIV_PIVSurf_Surface(1));
index_L=find(XPIV_LFV==0);
index_R=find(XPIV_LFV==length(PIV_match_Surface)+1);
index_RR=find(XPIV_LFV==XPIV_PIVSurf_Surface(end));

a = LFV_Surface_PIVMatched(1:index_LL-1);
b_L=LFV_Surface_PIVMatched(index_LL-1:index_L);
b = LFV_Surface_PIVMatched(index_LL:index_RR);
b_R=LFV_Surface_PIVMatched(index_R:index_RR);
c = LFV_Surface_PIVMatched(index_RR+1:end);
PS_Mod = PIVSurf_Surface;
PS_Mod_L=PS_Mod(1:length(b_L));
PS_Mod_R=PS_Mod(length(PS_Mod)-length(b_R)+1:end);

a=a+(mean(PS_Mod_L-b_L)); b_L=b_L+(mean(PS_Mod_L-b_L));
c=c+(mean(PS_Mod_R-b_R)); b_R=b_R+(mean(PS_Mod_R-b_R));

TZ_L = length(b_L);
TZ_R = length(b_R);
PS_Mod(1:TZ_L) = 1/(TZ_L-1) * ((0:TZ_L-1) .* PIVSurf_Surface(1:TZ_L) + (TZ_L-1:-1:0) .* b_L(1:TZ_L));
PS_Mod(end-TZ_R+1:end) = 1/(TZ_R-1) * ((0:TZ_R-1) .* b_R(end-TZ_R+1:end) + (TZ_R-1:-1:0) .* PIVSurf_Surface(end-TZ_R+1:end));
        
LFV_Surface_Combo_Surface = [a PS_Mod c];
LFV_Surface_Combo_Surface(length(a)-350:length(a)+150) = smooth(LFV_Surface_Combo_Surface(length(a)-350:length(a)+150),0.12);
LFV_Surface_Combo_Surface(end-length(c)+1-500:end-length(c)+50) = smooth(LFV_Surface_Combo_Surface(end-length(c)+1-500:end-length(c)+50),0.12);
% i=isnan(LFV_Surface_Combo_Surface);
% LFV_Surface_Combo_Surface=despike_fab(LFV_Surface_Combo_Surface);
% LFV_Surface_Combo_Surface(i)=NaN;
% LFV_Surface_Combo_Surface(index_L+1:index_R-1)=PIV_match_Surface;

% i=~isnan(LFV_Surface_Combo_Surface);
% LFV_Surface_Combo_Surface=LFV_Surface_Combo_Surface(i);
% XPIV_LFV=XPIV_LFV(i);
% LFV_Surface_PIVMatched=LFV_Surface_PIVMatched(i);

BadFrameLFVSurf = imSurf_LFV.badFrameBool;
if(max(LFV_Surface_Combo_Surface) > size(PIV1,1) || min(LFV_Surface_Combo_Surface) < 1)
    BadFrameLFVSurf=1;
end




        
        