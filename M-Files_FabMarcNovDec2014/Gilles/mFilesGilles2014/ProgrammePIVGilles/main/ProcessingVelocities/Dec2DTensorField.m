clear all
close all
%% Paramètres d'entrée
% exp_name = 'LC2_dt15ms_cc_1';
Res = 6.4101e-05; %m/pixel A VERIFIER
deltaT = 7d-3; %s
GrdSpc = [128 64 32 16 8]; 
start_image_number = 363;
end_image_number = 363;

%% Paramètres intrinsèques
% path = '\\beo\data\';
num_of_digits = 4;
% save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];

for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    %     load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt']);
    load(['E:\ComputedVelocities\20140804\LC1_1\ExpLC1_1_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
    
%     U=compVel.delx_ints*Res/deltaT;
%     V=compVel.dely_ints*Res/deltaT;
    U=compVel1.delta_x*Res/deltaT;
    V=compVel1.delta_z*Res/deltaT;
    %% Calcul des differentes composantes du tenseur du gradient des vitesses
    for i = 2:size(V,1)-1 
        for j = 2:size(V,2)-1
            a(i,j)=(U(i,j+1)-U(i,j-1))/(2*GrdSpc(end)*Res);
            b(i,j)=(U(i+1,j)-U(i-1,j))/(2*GrdSpc(end)*Res);
            c(i,j)=(V(i,j+1)-V(i,j-1))/(2*GrdSpc(end)*Res);
            d(i,j)=(V(i+1,j)-V(i-1,j))/(2*GrdSpc(end)*Res);
        end
    end
    %% Decomposition du tenseur en des termes de derfmation isotropique pure, rotation pure et rotation anisotrope.
    DefIso = (a+d)/2;
    Rotation = (c-b)/2;
    CisailCoef = sqrt((a-d).^2+(b+c).^2)/2;
    CisailRot = atan((b+c)./(a-d));
    
end


figure, imagesc(compVel1.delta_z.*compVel1.mask)
hold on, quiver(compVel1.delta_x.*compVel1.mask,compVel1.delta_z.*compVel1.mask)
