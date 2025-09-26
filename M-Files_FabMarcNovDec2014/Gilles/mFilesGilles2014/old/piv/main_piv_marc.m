% Main main_piv_marc
%% Object:
% computes velocities in the air above the water, and saves results to file
%% Functions called: 
% 1. fuseImages
% 2. findSurface
% 3. computeVelocities
%% Result: 
% experiment name, image pair number, delta_x (horizontal displacements),
% delta_z (vertical displacements), correlations, last level reached, mask
% used
%% Author:
% Marc Buckley
%% Last update:
% 05/09/2013
%%
tic,
clear all
close all
%% function arguments
exp_name = '4';
path = 'D:\Gilles_LC_PIV\testImages\RawImages\';
num_of_digits = 4;
piv_res = 47.4; %um/pix
vec_res = 4 * piv_res;
pivsurf_res = 99.6; %um/pix
piv_delta_t = 1/7.2; %sec
%% main arguments
save_path = 'D:\Gilles_LC_PIV\testImages\RawImages\';
start_image_pair_number = 2247;
end_image_pair_number = 2247;  %!!!!!!!!!!!!!!!!!!
%
for image_pair_number = start_image_pair_number:end_image_pair_number
%     image_pair_number = start_image_pair_number%:end_image_pair_number
    disp(['pair ' num2str(image_pair_number) ' finding surface...']);
    %% fuse first frames of pair
    image_letter = 'a';
    fusedIm = fuseImages(exp_name, [path 'Piv1\'],[path 'Piv2\'], image_pair_number, num_of_digits, image_letter);
    IM1 = fusedIm.fused_im;
    toc % !!!!!!!!!!!!!!
    %% find surface on first fused frame and generate mask for PIV computation
    imSurf = findSurface(exp_name, [path 'Pivsurf\'], image_pair_number, num_of_digits, image_letter, IM1, piv_res, pivsurf_res);
    mask = imSurf.mask;
    x_s = imSurf.x_s;
    z_s = imSurf.z_s;
    toc % !!!!!!!!!!!!!!
    disp(['pair ' num2str(image_pair_number) ' surface found, now computing velocities...']);
    %% fuse second frames of pair
    image_letter = 'b';
    fusedIm = fuseImages(exp_name, [path 'Piv1\'],[path 'Piv2\'], image_pair_number, num_of_digits, image_letter);
    IM2 = fusedIm.fused_im;
    %% compute velocities
    compVel =  computeVelocities(exp_name, image_pair_number, IM1, IM2, mask); 
    %% save to file
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel', 'x_s', 'z_s');
    disp(['pair ' num2str(image_pair_number) ' done.']);
end
toc