function [PIV3_CamAngle_Corrected] = Correct_Angle_PIV(PIV3)
% Correct PIV Surf Camera Angle

%% Correct residual distortion through polynomial
% Load imagePoints from PIVSURF_match (see PIVSURF_match.m in PIV
% calibration)
load PIVSURF_match.mat U1 X1

%%% This polynomial transformation also corrects for residual errors from camera
%%% distortion
T1 = fitgeotrans(U1,X1,'polynomial',4);
PIV3(isnan(PIV3)) = 0;
PIV3_CamAngle_Corrected = imwarp(PIV3,T1,'linear');

%% Further correction to get straight plumblines and flat surface through projective
PL_L(1,:) = [160,5336];
PL_L(2,:) = [4350,5354];
%cc_1 = polyfit(PL_L(:,1),PL_L(:,2),1);
PL_R(1,:) = [160,830];
PL_R(2,:) = [4781,824];
%cc_2 = polyfit(PL_R(:,1),PL_R(:,2),1);
FS(1,:) = [5864,5360.5];
FS(2,:) = [5860,822.5];
%cc_3 = polyfit(FS(:,1),FS(:,2),1);
%figure;imagesc(PIV2_CamAngle_Corrected);hold on;plot(cc_1(1)*[1:7920,1:7920]+cc_1(2),'r',cc_2(1)*[1:7920,1:7920]+cc_2(2),'r',cc_3(1)*[1:7920,1:7920]+cc_3(2),'r')

%%% Points retrieved from plumblines and flat surface
UU1 = [PL_R(1,:) ; PL_L(1,:) ; FS];
XX1 = [ 160 825.25 ; 160 5348.25 ; 5862 5348.25 ; 5862 825.25 ];

%%% Projective transformation
%T1 = maketform('projective',U1,X1); % Creates spatial transformation struct
TT1 = fitgeotrans(UU1,XX1,'projective');
% for a two-dimensional projective transformation that map each row of U to
% the corresponding  row of X. The U and X arguments are each 4-by-2 matrix
% and define the corners of input and output quadrilaterals.  Note that  no
% three corners can be collinear.
PIV3_CamAngle_Corrected = imwarp(PIV3_CamAngle_Corrected,TT1,'linear');

% The  corrected PIV  surface image for
% the camera angle, which is saved as "CamAngle" in the workspace.
% % % PIV2_CamAngle_Corrected=PIV2_CamAngle_Corrected(1:end-15,:);

%%% Crop NaN values
UCrop = [7997 6046 ; 7959 102 ; 161 30 ; 53 6025]; % vertices of the valid values
XCrop = [7959 6025 ; 7959 102 ; 161 102 ; 161 6025]; % square of the valid values
PIV3_CamAngle_Corrected = imcrop(PIV3_CamAngle_Corrected, [min(XCrop(:,1)) min(XCrop(:,2)) XCrop(1,1)-XCrop(3,1) XCrop(1,2)-XCrop(2,2)]);

%%% Flip image to have wave from left to right
PIV3_CamAngle_Corrected = PIV3_CamAngle_Corrected';

end