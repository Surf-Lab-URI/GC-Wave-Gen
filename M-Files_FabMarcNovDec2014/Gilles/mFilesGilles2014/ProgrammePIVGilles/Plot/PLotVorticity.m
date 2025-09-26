tic,
clear all
close all

%% Parametres
exp_name = 'LC1_1';
ResPivReal = 6.4101e-05; % m/pixel
ResPiv=1;
DeltaT= 7e-3;
Sc = ResPivReal/DeltaT;
hResLfv = 4.5663e-4/ResPivReal; 
vResLfv = 4.4956e-4/ResPivReal; 
num_of_digits = 4;
IntrWndw = [128 64 32 16 8];
GrdSpc = [64 32 16 8 4];
n = 2048; m = 2048;
path = ['E:\ComputedVelocities\20140804\'];
start_image_pair_number = 300;
end_image_pair_number = 500;
% %  Remarque : dt=12.5ms
u = [173 1756 968]'; v = [765 763 1819]'; x = [1505 586 1041]'; y = [1080 1091 510]'; 
y=y-16;
tform = maketform('affine',[x y],[u v]);
clear u v x y

for image_pair_number = start_image_pair_number:end_image_pair_number
    load([path exp_name '\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
    u = compVel1.delta_x*Sc;
    v = compVel1.delta_z*Sc;
    Vorticity = curl(u,v);
    Vorticity_rs = interpTransfoInverseNew(Vorticity,SUMASK);
    Vorticity_rs_m = nanmean(Vorticity_rs,2);
    Vorticity_rs_turb = Vorticity_rs - repmat(Vorticity_rs_m,1,size(u,2));
    Vorticity_rs_turb_m2 = nanmean(Vorticity_rs_turb.^2,2);
    Vorticitym(:,image_pair_number-start_image_pair_number+1) = Vorticity_rs_m;
    VorticityTurbm2(:,image_pair_number-start_image_pair_number+1) = Vorticity_rs_turb_m2;
end

Vorticity.v=Vorticitym;
Vorticity.VorticityTurbm2 = VorticityTurbm2;
Vorticity.start_image_pair_number = start_image_pair_number;
Vorticity.end_image_pair_number = end_image_pair_number;

filename = ['\Exp' exp_name '_Vorticity'];
outfile = [path exp_name filename];
save(outfile, 'Vorticity');




%% Plot 
inter=[1:150];
xx=[Vorticity.start_image_pair_number:Vorticity.end_image_pair_number]/7.2-21;
fc=328/269; % rectifie la modification d'echelle due a la transformation
yy=inter*ResPivReal*GrdSpc(end)*fc;
figure, imagesc(xx,yy,Vorticity.v(inter,:).^2), colorbar

figure, imagesc(xx,yy,Vorticity.VorticityTurbm2(inter,:)), colorbar
