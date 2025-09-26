%% Evolution de l'énergie en fonction du nombre d'onde et du temps
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC3_dt25ms_1';
Res = 1.2700e-04; %m/pixel A VERIFIER
deltaT = 25d-3; %s
start_image_number = 320;
end_image_number = 450;
hauteur_surf_repos = 150; % Paramètre a ajuster, on considère que la surface est encore suffisamment plane aux instants qui nous intéressent.
largeure_bande = 15; % Largeur de la bande sur laquelle est fait le moyennage.
%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% Calcul du spectre spatial au temps i
for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    V=sqrt((compVel.delx_ints*Res/deltaT).^2+(compVel.dely_ints*Res/deltaT).^2); %Ici le calcul est fait sur la valeur absolue de la vitesse, ce choix est-il le plus judicieux?
    for i= hauteur_surf_repos: hauteur_surf_repos + largeure_bande
        M(i,:)=pwelch(V(i,:));
    end
    spectre.sp(:,image_pair_number-start_image_number+1)=mean(M,1);
end
spectre.start_image_number = start_image_number;
spectre.end_image_number = end_image_number;
spectre.hauteur_surf_repos = hauteur_surf_repos;
spectre.largeure_bande = largeure_bande;

filename = ['Exp' exp_name '_Spectre' ];
outfile = [save_path filename];
save(outfile, 'spectre')

%% Analyse
% figure, loglog(N(:,1))
% hold on, loglog(N(:,10),'r'), loglog(N(:,30),'g')