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
exp_name = 'LC1_dt25ms_1';
path = 'F:\data\ExpLC1_dt25ms_1\rawImages';
num_of_digits = 4;
piv_res = 40.0; %um/pix
vec_res = 4 * piv_res;
pivsurf_res = 120.0; %um/pix
piv_delta_t = 1/7.2; %sec

%% Angle and size cam correction
u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
%% main arguments
save_path = 'D:\Gilles_LC_PIV\testImages\RawImages\';
start_image_pair_number = 700;
end_image_pair_number = 700;  
%
for image_pair_number = start_image_pair_number:end_image_pair_number
%     image_pair_number = start_image_pair_number%:end_image_pair_number
    disp(['pair ' num2str(image_pair_number) ' finding surface...']);
    %% fuse first frames of pair
    image_letter = 'a';
    fusedImLC = fuseImagesLC(exp_name, [path '\Piv1\'],[path '\Piv2\'], image_pair_number, num_of_digits, image_letter,tform);
    IM1 = fusedImLC.fused_im;
    %% find surface on first fused frame and generate mask for PIV computation
    imSurfLC = findSurfaceLC(exp_name, [path '\Pivsurf\'], image_pair_number, num_of_digits, image_letter, IM1,tformSurf);
    mask = imSurfLC.mask;
    x_s = imSurfLC.x_s;
    z_s = imSurfLC.z_s;
    disp(['pair ' num2str(image_pair_number) ' surface found, now computing velocities...']);
    %% fuse second frames of pair
    image_letter = 'b';
    fusedImLC = fuseImagesLC(exp_name, [path '\Piv1\'],[path '\Piv2\'], image_pair_number, num_of_digits, image_letter, tform);
    IM2 = fusedImLC.fused_im;
    %% compute velocities
    compVel =  computeVelocities(exp_name, image_pair_number, IM1, IM2, mask); 
    %% save to file
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel', 'x_s', 'z_s');
    disp(['pair ' num2str(image_pair_number) ' done.']);
end
toc