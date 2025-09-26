tic,
clear all
close all


%% Parameter
exp_name = 'LC2_dt25ms_cc_1';
deltaT = 25d-3; %s

%%
path = '\\beo\data\';
num_of_digits = 4;

IntrWndw = [128 64 32 16 8]; % IntrWndw=[128 64 32 16 8];
GrdSpc = [64 32 16 8 4]; 

save_path = ['\\beo\data\Exp' exp_name '\ComputedVelocities\'];
start_image_pair_number = 801;
end_image_pair_number = 900;


u = [58 1839 28 1878]'; v = [735 735 1808 1808]'; x = [58 1839 58 1839]'; y =[735 735 1808 1808]'; 
tformPivcc = maketform('projective',[u v],[x y]);

u = [253 1642 283 1618]'; v = [1060 1060 1882 1882]'; x = [58 1839 58 1839]'; x=x-62; y =[735 735 1808 1808]';
tformPivsurfcc = maketform('projective',[u v],[x y]);

clear u v x y
%%
 for image_pair_number = start_image_pair_number:end_image_pair_number
% image_pair_number = 327

    %% Images aux temps 'a'
    image_letter='a';
     %% Correction pivcc et pivsurfcc
    pivcc = load([path 'Exp' exp_name '\RawImages\Pivcc\' 'Exp' exp_name '_Pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    pivsurfcc = load([path 'Exp' exp_name '\RawImages\Pivsurfcc\' 'Exp' exp_name '_Pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    pivcc = pivcc.imgPivcc;
    pivsurfcc = pivsurfcc.imgPivsurfcc;
    %% Pivcc
    imgPivcc_rect=imtransform(pivcc,tformPivcc,'XYScale',1);
    imgPivcc_rect_cut=imgPivcc_rect(1:2086,62:2015);
    
    IM1 = imgPivcc_rect_cut;
    
    %% Pivsurfcc
    imgPivsurfcc_rect=imtransform(pivsurfcc,tformPivsurfcc,'XYScale',1);
    imgPivsurfcc_rect_cut=imgPivsurfcc_rect(1265-783:end,369:1953+369);
      
    surf1=findSurface_lc(imgPivsurfcc_rect_cut);
    mask1 = nan(size(imgPivsurfcc_rect_cut));
    for i=1:size(imgPivsurfcc_rect_cut,2)
        mask1(round(surf1.z_s(i)):end,i) = 1;
    end
%     clear pivsurf_struc imgPivsurf imgPivsurf_t_cut
    
    %% Images aux temps 'b'
    image_letter='b';
     %% Correction pivcc et pivsurfcc
    pivcc = load([path 'Exp' exp_name '\RawImages\Pivcc\' 'Exp' exp_name '_Pivcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    pivsurfcc = load([path 'Exp' exp_name '\RawImages\Pivsurfcc\' 'Exp' exp_name '_Pivsurfcc_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
    pivcc = pivcc.imgPivcc;
    pivsurfcc = pivsurfcc.imgPivsurfcc;
    %% Pivcc
    imgPivcc_rect=imtransform(pivcc,tformPivcc,'XYScale',1);
    imgPivcc_rect_cut=imgPivcc_rect(1:2086,62:2015);
    
    IM2 = imgPivcc_rect_cut;
    
    %% Pivsurfcc
    imgPivsurfcc_rect=imtransform(pivsurfcc,tformPivsurfcc,'XYScale',1);
    imgPivsurfcc_rect_cut=imgPivsurfcc_rect(1265-783:end,369:1953+369);
      
    surf2=findSurface_lc(imgPivsurfcc_rect_cut);
    mask2 = nan(size(imgPivsurfcc_rect_cut));
    for i=1:size(imgPivsurfcc_rect_cut,2)
        mask2(round(surf2.z_s(i)):end,i) = 1;
    end
%     clear pivsurf_struc imgPivsurf imgPivsurf_t_cut

    %% piv calc
    compVel = PIV_FAB6_LfvOrb1_noOrb (IM1,IM2, mask1, mask2, IntrWndw, GrdSpc);
    
    filename = ['Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'compVel');
    disp(['pair ' num2str(image_pair_number) ' done.']);
    
 end
toc