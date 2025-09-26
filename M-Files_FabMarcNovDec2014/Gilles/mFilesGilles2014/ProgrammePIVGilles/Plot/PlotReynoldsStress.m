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
    ut = interpTransfoInverseNew(compVel1.delta_x*Sc,SUMASK);
    vt = interpTransfoInverseNew(compVel1.delta_z*Sc,SUMASK);
    U = nanmean(ut,2);
    V = nanmean(vt,2);
    u = ut-repmat(U,1,size(ut,2));
    v = vt-repmat(V,1,size(vt,2));
    
    u2m(:,image_pair_number-start_image_pair_number+1) = nanmean(u.^2,2);
    v2m(:,image_pair_number-start_image_pair_number+1) = nanmean(v.^2,2);
    uvm(:,image_pair_number-start_image_pair_number+1) = nanmean(u.*v,2);  
end    

RS.u2m = u2m;
RS.v2m = v2m;
RS.uvm = uvm;
RS.start_image_pair_number = start_image_pair_number;
RS.end_image_pair_number = end_image_pair_number;
filename = ['\Exp' exp_name '_RS'];
outfile = [path exp_name filename];
save(outfile, 'RS');
%% Plot 
inter=[1:150];
fc=328/269; % rectifie la modification d'echelle due a la transformation
xx=[RS.start_image_pair_number:RS.end_image_pair_number]/7.2-21;
yy=inter*ResPivReal*GrdSpc(end)*fc;
figure, imagesc(xx,yy,u2m(inter,:)), colorbar, caxis([0 0.005])
figure, imagesc(xx,yy,v2m(inter,:)), colorbar, caxis([0 0.005])
figure, imagesc(xx,yy,uvm(inter,:)), colorbar, caxis([-0.00075 0.00075])