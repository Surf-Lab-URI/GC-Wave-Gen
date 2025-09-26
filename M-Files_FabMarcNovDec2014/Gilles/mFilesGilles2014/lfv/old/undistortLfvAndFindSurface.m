% undistort lfv and find surface
%% Object:

%% Functions called:
%% Result:
%
%% Author:
% Marc Buckley
%% Last update:
% 07/02/2013
%%
tic,
clear all
close all
load('\\beo\mFiles\lfv\Calib_Results.mat')
% load('E:\data\PIV\Exp4\RawImages\Lfv\Exp4_Lfv_0008.mat');

%% user input
exp_name = '3';
start_image_pair_number = 0;
end_image_pair_number = 3023;

%%
path = ['E:\data\Exp' exp_name '\' 'RawImages\'];
num_of_digits = 4;
lfvWidth = 0.5598; %m
lfvHeight = 0.5125; %m
lfvWres = lfvWidth/2048; %m/pix
lfvHres = lfvHeight/2048;
lfvSWL =  2048 - 1.0125e+003 + 1; %pixel
%% main arguments
save_path = ['E:\data\Exp' exp_name '\' 'ProcessedLfv\'];
% start_image_pair_number = 0;
% end_image_pair_number = 5;
lfvSurfElevationWrtSWL = nan(end_image_pair_number+1,2048);
% badIm = char(end_image_pair_number+1,num_of_digits);
badImNumber = 0;
%
for image_pair_number = start_image_pair_number:end_image_pair_number
    % image_pair_number = start_image_pair_number%:end_image_pair_number
    
    %%
    
    %% UNDISTORT THE IMAGE:
    disp(['pair ' num2str(image_pair_number) ' undistorting image...']);
    load([path 'Lfv\Exp' exp_name '_Lfv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
    I = imgLfv;
    KK = [fc(1) alpha_c*fc(1) cc(1) ; 0 fc(2) cc(2) ; 0 0 1];
    
    %%% Compute the new KK matrix to fit as much data in the image (in order to
    %%% accomodate large distortions:
    r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
    dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
    fc_new = dist_amount * fc;
    
    KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];
    
    [I2] = rect(I,eye(3),fc,cc,kc,alpha_c,KK_new);
    
    %% find surface
    disp(['pair ' num2str(image_pair_number) ' finding surface...']);
%     imSurf1 = findSurfaceLfvFabMarc(I2, thresh);
    imSurf1 = findSurfaceLfv(I2);
    
    
    z_s = (2048 - imSurf1.z_s + 1 - lfvSWL) * lfvHres;
    
    lfvSurfElevationWrtSWL(image_pair_number+1,:) = z_s;
    
%     if imSurf1.badImColCounter>50
%         badImNumber = badImNumber + 1;
%         disp(['pair ' num2str(image_pair_number) ' is bad.']);
%         badIm(badImNumber,:) = sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number);
%         figure('name', sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number),'NumberTitle','off'), imagesc(imgLfv), colormap(bone), hold on, plot(imSurf1.x_s, imSurf1.z_s, 'r')
% %         figure('name', sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number),'NumberTitle','off'), imagesc(imgLfv), colormap(bone), hold on, plot(imSurf1.x_s, imSurf1.z_s_raw, 'r')
%         keyboard
%     end
    

% figure('name', sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number),'NumberTitle','off'), imagesc(imgLfv), colormap(bone), hold on, plot(imSurf1.x_s, imSurf1.z_s, 'r')
% figure('name', sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number),'NumberTitle','off'), imagesc(imgLfv), colormap(bone), hold on, plot(imSurf1.x_s, imSurf1.z_s_raw, 'r')
% keyboard
    disp(['pair ' num2str(image_pair_number) ' done.']);
end
x_s = (imSurf1.x_s-1) * lfvWres;
t_s = 0:1/7.2:end_image_pair_number * 1/7.2;

filename = ['Exp' exp_name '_lfvSurfElevationWrtSWL'];
outfile = [save_path filename];
%%
% badProfiles = badIm;
% badProfiles = [6 20 26 222 227 236 242 443 452 654 668 674 884 890 1090 1099 1104];
%%
save(outfile, 't_s','x_s', 'lfvSurfElevationWrtSWL');
toc

