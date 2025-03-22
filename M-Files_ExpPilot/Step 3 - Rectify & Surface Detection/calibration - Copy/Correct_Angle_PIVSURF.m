function [PIV2_CamAngle_Corrected] = Correct_Angle_PIVSURF(PIV2)

%%%%%%%% Correct PIV Surf Camera Angle %%%%%%%%

%% Correct residual distortion through polynomial
% Load imagePoints from PIVSURF_match (see PIVSURF_match.m in PIV
% calibration)
load PIV_match.mat U1 X1 Xt

%%% This polynomial transformation also corrects for residual errors from camera
%%% distortion
T1 = fitgeotrans(U1,X1,'polynomial',4);
PIV2(isnan(PIV2)) = 0;
PIV2_CamAngle_Corrected2 = imwarp(PIV2,T1,'linear');
PIV2_CamAngle_Corrected2 = imtranslate(PIV2_CamAngle_Corrected2, [(X1(1,1)-Xt(1,1)) (X1(1,2)-Xt(1,2))] );

%% Correction to get straight plumblines and flat surface through projective
% % % % Without polynomial correction
% % % % % % CamAngle = PIV2;
% % % % % % PL_L = [123 8 ; 181.5 875];
% % % % % % %cc_1 = polyfit(PL_L(:,1),PL_L(:,2),1);
% % % % % % PL_R = [3945.5 914 ; 4008 8];
% % % % % % %cc_2 = polyfit(PL_R(:,1),PL_R(:,2),1);
% % % % % % FS = [224 1508.25 ; 3904.5 1507.5];
% % % % % % %cc_3 = polyfit(FS(:,1),FS(:,2),1);
% % % 
% % % % With polynomial correction
% % % CamAngle = PIV2_CamAngle_Corrected2;
% % % PL_L = [373 5 ; 368 849];
% % % %cc_1 = polyfit(PL_L(:,1),PL_L(:,2),1);
% % % PL_R = [4087 888 ; 4082 7];
% % % %cc_2 = polyfit(PL_R(:,1),PL_R(:,2),1);
% % % FS = [364 1501.9 ; 4090.5 1501.8];
% % % %cc_3 = polyfit(FS(:,1),FS(:,2),1);
% % % 
% % % %figure;imagesc(PIV2);hold on;plot(1:4096,cc_1(1)*[1:4096]+cc_1(2),'r',1:4096,cc_2(1)*[1:4096]+cc_2(2),'r',1:4096,cc_3(1)*[1:4096]+cc_3(2),'r')
% % % 
% % % %%% Points retrieved from plumblines and flat surface
% % % UU1 = [PL_L(1,:) ; FS ; PL_R(2,:)];
% % % XX1 = [ 368.5 6 ; 368.5 1501.9 ; 4086 1501.9 ; 4086 6 ];
% % % 
% % % %%% Projective transformation to correct PIVSURF Camera Angle
% % % %T1 = maketform('projective',U1,X1); % Creates spatial transformation struct
% % % TT1 = fitgeotrans(UU1,XX1,'projective');
% % % % for a two-dimensional projective transformation that map each row of U to
% % % % the corresponding  row of X. The U and X arguments are each 4-by-2 matrix
% % % % and define the corners of input and output quadrilaterals.  Note that  no
% % % % three corners can be collinear.
% % % PIV2_CamAngle_Corrected = imwarp(CamAngle,TT1,'linear'); %'cubic' instead of 'linear'?

% Flip image to have wave from left to right
PIV2_CamAngle_Corrected = fliplr(PIV2_CamAngle_Corrected2);

end