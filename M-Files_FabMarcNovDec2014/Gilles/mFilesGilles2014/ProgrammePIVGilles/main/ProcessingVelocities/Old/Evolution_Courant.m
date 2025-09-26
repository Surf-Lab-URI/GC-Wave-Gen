%% Evolution de du courant suivant l'axe x en fonction de la profondeur (definie a partir de SU) et du temps
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC3_dt7ms_3';
Res = 1.2700e-04; %m/pixel A VERIFIER
deltaT = 7d-3; %s
start_image_number = 355;
end_image_number = 365;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% Calcul du spectre spatial au temps i
for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    compVel.courant=nanmean(compVel.delx_ints.*compVel.MASK,2);
    courant.cour(:,image_pair_number-start_image_number+1)=compVel.courant;
end
courant.start_image_number = start_image_number;
courant.end_image_number = end_image_number;

filename = ['Exp' exp_name '_Courant' ];
outfile = [save_path filename];
save(outfile, 'courant')

%% Analyse
% for i=1:size(courant.cour,1)
%     smoothnCour(i,:)=smoothn(courant.cour(i,:),10);
% end
% 
% 
% % Video
% for image_pair_number=start_image_number:end_image_number
%     plot(courant.cour(:,image_pair_number-start_image_number+1))
%     G(image_pair_number)=getframe;
% end
% movie(G,3,1)