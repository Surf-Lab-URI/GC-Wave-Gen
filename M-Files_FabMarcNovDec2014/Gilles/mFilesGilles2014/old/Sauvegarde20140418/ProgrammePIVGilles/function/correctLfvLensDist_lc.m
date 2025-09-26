function imgLfv_LD = correctLfvLensDist(imgLfv)


%% correct image for lens distortion
load('\\beo\mFiles\lfv\Calib_Results.mat')
I = imgLfv;

% KK = [fc(1) alpha_c*fc(1) cc(1) ; 0 fc(2) cc(2) ; 0 0 1];
    
%%% Compute the new KK matrix to fit as much data in the image (in order to
%%% accomodate large distortions:
% r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
fc_new = dist_amount * fc;

KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

[I2] = rect(I,eye(3),fc,cc,kc,alpha_c,KK_new);

imgLfv_LD = I2;