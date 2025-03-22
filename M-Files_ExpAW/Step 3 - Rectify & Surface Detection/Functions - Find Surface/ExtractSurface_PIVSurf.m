function [XPos,YPos,PIVSurf_PIVMatch,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIV_match_Surface,BadFramePIVSurf] = ExtractSurface_PIVSurf(PIVSURF_CR1,PIV1,expName)

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

%% Fabrice's old PIVSurf (for URI2018 experiments)
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

%%% Procedure for expts 1-48 (done in June '22)
if str2double(expName)<49
    % Polynomial transformation to correct for residual errors from camera
    % distortion
    load PIV_match.mat U1 X1 Xt
    T1 = fitgeotrans(U1,X1,'polynomial',4);
    PIVSURF_CamAngle = PIVSURF_CR1;
    PIVSURF_CamAngle(isnan(PIVSURF_CamAngle)) = 0;
    PIVSURF_CamAngle = imwarp(PIVSURF_CamAngle,T1,'linear');
    PIVSURF_CamAngle = imtranslate(PIVSURF_CamAngle, [(X1(1,1)-Xt(1,1)) (X1(1,2)-Xt(1,2))] );
    PIVSURF_CamAngle = fliplr(PIVSURF_CamAngle);
    
    % Find transformation to match pixel's resolution on PIV image
    P = load('PIV_PIVsurf_matching_points.mat','Xp'); %,'Up','PIV_iP','PIVSURF_iP','PIVSURF_iP2');
    Up2 = [-X1(:,1)+size(PIVSURF_CamAngle,2)+1,X1(:,2)];
    Up2([1:26,end-8*26+1:end],:) = [];
    Up2 = flipud(reshape(permute(reshape(Up2,[26,31,2]),[2,1,3]),26*31,2));
    Up = [ Up2(end-30,:) ; Up2(end,:) ; Up2(31,:) ; Up2(1,:) ]; % The coordinate
    % of the points in the PIVSurf image (the grid calibration image).
    Xp = P.Xp; % The coordinates of the points on PIV image. The physical location  of these points are
    % the same as PIVSurf images (Up). The calibration grid was used to find out
    % the exact same locations for these two images.
    T2_bis = maketform('projective',Up,Xp);
    
    % PIVSurf image at PIV resolution --> this is BEFORE further rotation and
    % translation!!!!!
    [Resized_PIVSurf_bis,XPos,YPos] =  imtransform(PIVSURF_CamAngle,T2_bis,'XYScale',1);
    
    % Cropped image to fit PIV size
    PIVSurf_PIVMatch_pre = Resized_PIVSurf_bis(1:round(size(PIV1,1)-YPos(1)),round(1-XPos(1)):round(size(PIV1,2)-XPos(1)));
    PIVSurf_PIVMatch_pre = [zeros(round(YPos(1)),size(PIVSurf_PIVMatch_pre,2));PIVSurf_PIVMatch_pre];
    
    %% Find surface on resized PIVSurf
    
    % Step 0: Transform Resized PIVSurf in order to match PIVSurf_PIVMatch
    % (i.e., construct image to rotate of Alpha around the same point and
    % translate according to DX2)
    Resized_PIVSurf_bis2 = [ones(round(YPos(1)),size(Resized_PIVSurf_bis,2));Resized_PIVSurf_bis];
    Cp = [size(PIVSurf_PIVMatch_pre,1)/2,size(PIVSurf_PIVMatch_pre,2)/2-round(XPos(1))]; % image center with PIVSurf_PIVMatch coordinates
    C = [size(Resized_PIVSurf_bis2,1)/2,size(Resized_PIVSurf_bis2,2)/2]; % image center with Resized_PIVSurf_bis2 coordinates
    DeltaCols = abs((Cp(2)-C(2))*2); % number of Columns to be added to make C coordinates equal to Cp
    DeltaRows = abs((Cp(1)-C(1))*2); % number of Rows to be added to make C coordinates equal to Cp
    R2 = [];
    if Cp(2)-C(2)>0 % if x(Cp) is bigger of x(C), add columns on the left
        R2 = [2*ones(DeltaCols,size(Resized_PIVSurf_bis2,2));Resized_PIVSurf_bis2];
    else % otherwise, add columns on the left
        R2 = [Resized_PIVSurf_bis2;2*ones(DeltaCols,size(Resized_PIVSurf_bis2,2))];
    end
    if Cp(1)-C(1)>0 % if z(Cp) is bigger of z(C), add rows on the top
        R2 = [2*ones(size(R2,1),DeltaRows),R2];
    else % otherwise, add rows on the bottom
        R2 = [R2,2*ones(size(R2,1),DeltaRows)];
    end
    % figure;imagesc(R2)
    
    %%%% To be adjusted for all tests (these are for Movie1 to 12)
    Alpha = rad2deg(9.15/5924); % further rotation checked with pre-test
    DX2 = 25.3894; % further translation after Alpha-rotation
    %%%%
    
    R2_bis = imrotate(R2,-Alpha,'nearest','crop');
    R2_bis(1:DeltaCols,:) = [];R2_bis(:,end-DeltaRows+1:end)=[]; % PIVSurf at
    
    % PIVSurf which matches PIV resolution and size
    % PIVSurf image matching PIV --> IT FITS PERFECTLY WITH PIV_match_Surface!!
    PIVSurf_PIVMatch = imtranslate(R2_bis(1:size(PIV1,1),round(1-XPos(1)):round(size(PIV1,2)-XPos(1))),[0,-DX2+1]);
    
    % Full-size PIVSurf which matches PIV resolution
    R2_bis2 = imtranslate(R2_bis(1:size(PIV1,1),:),[0,-DX2+1]); % full-size PIVSurf which matches PIV coordinates
end

%%% Procedure for expts 49-72 (done in December '22)
if str2double(expName)>48
    
%     PIVSURF_CamAngle = imrotate(PIVSURF_CR1,-rad2deg(4/4000));
    PIVSURF_CamAngle = PIVSURF_CR1;
    
    % Find transformation to match pixel's resolution on PIV image
    Up = [ 991 1931 ; 3177 1932 ; 3256 225 ; 916 227 ]  ; % The coordinate
    % of the points in the PIVSurf image (the grid calibration image).
    Xp = [ 5265 6708 ; 241 6712 ; 234 2861 ; 5280 2860 ]; % The coordinates of the points on PIV image. The physical location  of these points are
    % the same as PIVSurf images (Up). The calibration grid was used to find out
    % the exact same locations for these two images.
    T2_bis = maketform('projective',Up,Xp);
    
    % PIVSurf image at PIV resolution --> this is BEFORE further rotation and
    % translation!!!!!
    [Resized_PIVSurf_bis,XPos,YPos] =  imtransform(PIVSURF_CamAngle,T2_bis,'XYScale',1);
    
    % Cropped image to fit PIV size
    PIVSurf_PIVMatch_pre = Resized_PIVSurf_bis(1:round(size(PIV1,1)-YPos(1)),round(1-XPos(1)):round(size(PIV1,2)-XPos(1)));
    PIVSurf_PIVMatch = [zeros(round(YPos(1)),size(PIVSurf_PIVMatch_pre,2));PIVSurf_PIVMatch_pre];
    
    % Step 0: rotation to have horizontal surface
    R2_bis2 = [zeros(round(YPos(1)),size(Resized_PIVSurf_bis,2));Resized_PIVSurf_bis(1:round(size(PIV1,1)-YPos(1)),:)];
end


% Step 1: Define filter S to enhance surface contrast
%%% S filter computed on PIVSurf_PIVmatch image
F1 = 4000; %5000;
F2 = 6500;
if str2double(expName)>24 && str2double(expName) < 37
    F2 = 7000;
end
if str2double(expName)>48
    F2 = 7000;
end
% % % Lint = 50; %500???       150? 250?
% % % if str2double(expName)>36
% % %     Lint = 250;
% % % end
% % % S_T = R2_bis2(F1+Lint+1:F2,:); %S_T = imfilter(S_T,fspecial('gaussian',64,16),'replicate');
% % % S_B = R2_bis2(F1+1:F2-Lint,:); %S_B = imfilter(S_B,fspecial('gaussian',64,16),'replicate');
% % % 
% % % SS = S_T-S_B;%SS(abs(SS)<10)=0;%figure;imagesc(SS)
% % % S_pos2 = SS;S_pos2(S_pos2<0)=0;
% % % S_neg2 = SS;S_neg2(S_neg2>0)=0;
% % % S = [zeros(F1+Lint,size(S_neg2,2));S_pos2+[S_neg2(Lint+1:end,:);zeros(Lint,size(S_neg2,2))];zeros(size(R2_bis2,1)-F2,size(S_neg2,2))];
% % % % S = [zeros(F1+Lint,size(S_neg,2));S_pos2;zeros(size(R2_bis2,1)-F2,size(S_neg,2))];
% % % S(S<0)=0;
% % % S_old = S;
% % % S_old(S_old>200)=100;

Lint = 500;
if str2double(expName) > 36
    Lint = 250;
end
S_T = R2_bis2(F1+Lint+1:F2,:); S_T = imfilter(S_T,fspecial('gaussian',64,16),'replicate');
S_B = R2_bis2(F1+1:F2-Lint,:); S_B = imfilter(S_B,fspecial('gaussian',64,16),'replicate');
SS = S_T-S_B;
S_pos2 = SS;S_pos2(S_pos2<0)=0;
S_neg2 = SS;S_neg2(S_neg2>0)=0;
S = [zeros(F1+Lint,size(S_neg2,2));S_pos2;zeros(size(R2_bis2,1)-F2,size(S_neg2,2))];
S(S>200)=80;
if str2double(expName) > 36
    for i = 1101:7000
        S(1:4250,i) = S(4251,i);
        S(6501:end,i) = S(6500,i);
    end
    for i = 7001:size(S,2)
        S(6501:end,i) = S(6500,i);
    end
    S(1:4250,:) = imfilter(S(1:4250,:),fspecial('gaussian',256,64),'replicate');
    S(6501:end,:) = imfilter(S(6501:end,:),fspecial('gaussian',256,64),'replicate');
end


% Step 1: find surface on PIVSurf image with original size
X = 1:size(PIVSURF_CamAngle,2);
Xgood = 1492:9310; % good points where image is still bright enough to detect surface 9300
if str2double(expName)>48 && str2double(expName)<55
    Xgood = 800:9310;
end
% Xgood = 1501:9310; % good points where image is still bright enough to detect surface 9300
X = Xgood;
% % % [imSurf_PIVSURF] = FindSurface(S_old,5,1,expName);
[imSurf_PIVSURF] = FindSurface(S,5,1,expName);
PIVSurf_Surface_Raw = despike_fab(imSurf_PIVSURF.surface(X));
PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
IntFilt = 100000; % stronger wind, lower filter
if ismember(str2double(expName),[1:3,37:39])
    IntFilt = 75000;
elseif ismember(str2double(expName),[4:6,13:15,16:18,25:27,28:30,40:42])
    IntFilt = 50000;
elseif ismember(str2double(expName),[7:9,19:21,31:33,43:45])
    IntFilt = 25000;
elseif ismember(str2double(expName),[10:12,22:24,34:36,46:48])
    IntFilt = 1000;
end
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, IntFilt); %1000 for expName ~= 1,13,25,37
% PIVSurf_Surface = SP.coefs; %%% use for flat surface
% PIVSurf_Surface = polyval(polyfit([X(1),X(end)],PIVSurf_Surface,1),X);  %%% use for flat surface
PIVSurf_Surface_2 = SP.coefs(2:end-1);

% X-coordinate in PIV resolution for detected surface
X_PIVSURF2 = X+XPos(1);
XPIV_PIVSurf_Surface = round(X_PIVSURF2(10)):round(X_PIVSURF2(end-10)); % :round(X_PIVSURF2(end)); % rounded coordinates 
% XPIV_PIVSurf_Surface = round(X_PIVSURF2(1)):round(X_PIVSURF2(end-10)); % :round(X_PIVSURF2(end)); % rounded coordinates 
PIVSurf_Surface = interp1(X_PIVSURF2,PIVSurf_Surface_2,XPIV_PIVSurf_Surface,'pchip','extrap'); % interpolation on rounded coordinates

%%% Further rotation to correct for different water depth
% Case 1: w.d. 8 cm ---> expt 1-12, no corrections

% Case 2: w.d. 6 cm ---> expt 13:24
if str2double(expName)>12 && str2double(expName)<25
    Alpha2 = -4.5/size(PIVSurf_PIVMatch,2); % further rotation
    DX2 = -3;
    MM = [cos(atan(Alpha2)) sin(atan(Alpha2)); -sin(atan(Alpha2)) cos(atan(Alpha2))];
    Rot = MM*[X_PIVSURF2;PIVSurf_Surface_2]; %MM*[XPIV_PIVSurf_Surface;PIVSurf_Surface];
    PIVSurf_Surface = interp1(Rot(1,:),Rot(2,:),XPIV_PIVSurf_Surface)-DX2; 
end

% Case 3: w.d. 6 cm ---> expt 25:36
if str2double(expName)>24 && str2double(expName)<37
    Alpha2 = -6/size(PIVSurf_PIVMatch,2); % further rotation
    DX2 = -5;
    MM = [cos(atan(Alpha2)) sin(atan(Alpha2)); -sin(atan(Alpha2)) cos(atan(Alpha2))];
    Rot = MM*[X_PIVSURF2;PIVSurf_Surface_2]; %MM*[XPIV_PIVSurf_Surface;PIVSurf_Surface];
    PIVSurf_Surface = interp1(Rot(1,:),Rot(2,:),XPIV_PIVSurf_Surface)-DX2;
end

% Case 4: w.d. ~70 cm ---> expt 37:48
if str2double(expName)>36 && str2double(expName)<49
    Alpha2 = 6/size(PIVSurf_PIVMatch,2); % further rotation
    DX2 = 5;
    MM = [cos(atan(Alpha2)) sin(atan(Alpha2)); -sin(atan(Alpha2)) cos(atan(Alpha2))];
    Rot = MM*[X_PIVSURF2;PIVSurf_Surface_2]; %MM*[XPIV_PIVSurf_Surface;PIVSurf_Surface];
    PIVSurf_Surface = interp1(Rot(1,:),Rot(2,:),XPIV_PIVSurf_Surface)-DX2;
end

% Case 5: w.d. 70.0 cm ---> expt 49:60
if str2double(expName)>48 && str2double(expName)<61
    Alpha2 = 7/length(PIVSurf_Surface); % further rotation
    DX2 = -6;
    MM = [cos(atan(Alpha2)) sin(atan(Alpha2)); -sin(atan(Alpha2)) cos(atan(Alpha2))];
    Rot = MM*[X_PIVSURF2;PIVSurf_Surface_2]; %MM*[XPIV_PIVSurf_Surface;PIVSurf_Surface];
    PIVSurf_Surface = interp1(Rot(1,:),Rot(2,:),XPIV_PIVSurf_Surface)-DX2;
end

% Case 6: w.d. ~68.9 cm ---> expt 61:72
if str2double(expName)>60 && str2double(expName)<73
    Alpha2 = 7/length(PIVSurf_Surface); % further rotation
    DX2 = -6;
    MM = [cos(atan(Alpha2)) sin(atan(Alpha2)); -sin(atan(Alpha2)) cos(atan(Alpha2))];
    Rot = MM*[X_PIVSURF2;PIVSurf_Surface_2]; %MM*[XPIV_PIVSurf_Surface;PIVSurf_Surface];
    PIVSurf_Surface = interp1(Rot(1,:),Rot(2,:),XPIV_PIVSurf_Surface)-DX2;
end

%%%%%%%%%%%%%%%

% PIVSurf surface matching PIV coordinates
IndX = [find(XPIV_PIVSurf_Surface==1),find(XPIV_PIVSurf_Surface==size(PIV1,2))];
PIV_match_Surface = PIVSurf_Surface(IndX(1):IndX(end));

%% Check bad frames
BadFramePIVSurf = imSurf_PIVSURF.badFrameBool;
if (max(PIV_match_Surface) > size(PIV1,1) || min(PIV_match_Surface) < 1)
    BadFramePIVSurf = 1;
end

