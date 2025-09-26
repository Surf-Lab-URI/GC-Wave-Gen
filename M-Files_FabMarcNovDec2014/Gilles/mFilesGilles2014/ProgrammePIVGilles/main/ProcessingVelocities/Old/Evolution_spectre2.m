%% Evolution de l'énergie en fonction du nombre d'onde et du temps
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC2_dt15ms_cc_1';
Res = 1.2700e-04; %m/pixel A VERIFIER
deltaT = 7d-3; %s
GrdSpc = [128 64 32 16 8]; 
start_image_number = 150;
end_image_number = 150;
hauteur_surf_repos = 250; % Paramètre a ajuster, on considère que la surface est encore suffisamment plane aux instants qui nous intéressent.
largeur_bande = 15; % Largeur de la bande sur laquelle est fait le moyennage.
%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% Calcul du spectre spatial au temps i
% for image_pair_number=start_image_number:end_image_number
%     disp(image_pair_number);
image_pair_number=start_image_number
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_medfilt']);
    V=compVel.delx_ints*Res/deltaT;
    for i= hauteur_surf_repos: hauteur_surf_repos + largeur_bande
        f = fft(V(i,:)); %Fourier Modal decomposition
        fa = 2*abs(f)/length(f); %amplitude of mode
        E(i,:)=fa;
        k = [0:(size(V,2)-1)]*pi/(((size(V,2)-1)*GrdSpc(end)*Res)/2);  
    end
    result=mean(E,1);
    spectre.sp(1:size(V,2)/2,image_pair_number-start_image_number+1)=result(1:size(V,2)/2);
% end
spectre.start_image_number = start_image_number;
spectre.end_image_number = end_image_number;
spectre.hauteur_surf_repos = hauteur_surf_repos;
spectre.largeur_bande = largeur_bande;

% filename = ['Exp' exp_name '_Spectre' ];
% outfile = [save_path filename];
% save(outfile, 'spectre')

%% Analyse
figure, loglog(k(1:size(V,2)/2),spectre.sp(:,1))
% hold on, loglog(spectre.sp(:,75),'r'), loglog(spectre.sp(:,150),'g')