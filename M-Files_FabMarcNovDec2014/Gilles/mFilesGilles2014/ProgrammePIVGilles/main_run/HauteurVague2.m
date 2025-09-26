tic,
clear all
close all

%% Parametres
exp_name = 'LC2_1';
ResPivReal = 6.4101e-05; % m/pixel
ResPiv=1;
hResLfv = 4.5663e-4/ResPivReal; 
vResLfv = 4.4956e-4/ResPivReal; 
num_of_digits = 4;
IntrWndw = [128 64 32 16 8];
GrdSpc = [64 32 16 8 4];
n = 2048; m = 2048;
path = ['E:\data\20140804\LC2_1\Movie4_Scene5'];
save_path = ['E:\ComputedVelocities\20140804\LC2_1\'];
start_image_pair_number = 1;
end_image_pair_number = 700;
% %  Remarque : dt=12.5ms
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]'; 
y=y-16;
tform = maketform('affine',[x y],[u v]);
clear u v x y

for image_pair_number = start_image_pair_number:end_image_pair_number
    image_pair_number
    rawFrame = [path '_lfv_' sprintf(['%0' num2str(num_of_digits-1) 'd'], image_pair_number) '.raw'];
    matFrame = saveToMatSingleFrame(rawFrame,n,m);
    Lfv= rot90(matFrame.img,2);
    Lfv_LD = correctLfvLensDist_lc(Lfv);
    Lfv_cut=Lfv_LD(401:1600,401:1600); % Seul une partie de l'image est expoitable
    surfLfv=findSurface_lc(Lfv_cut);
    HVague(image_pair_number-start_image_pair_number+1)= max(surfLfv.z_s)- min(surfLfv.z_s);
end

filename = ['Exp' exp_name '_HVague'];
outfile = [save_path filename];
save(outfile, 'HVague');