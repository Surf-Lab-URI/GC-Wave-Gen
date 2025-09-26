%% Evolution de l'énergie cinetique en fonction de la profondeur (definie a partir de SU) et du temps
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC2_dt7ms_3';
Res = 40.00e-06; %m/pixel A VERIFIER
deltaT = 7d-3; %s
start_image_number = 250;
end_image_number = 650;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% Calcul du spectre spatial au temps i
for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    energCin.ec(:,image_pair_number-start_image_number+1)=1/2*(compVel.courantx.^2+compVel.couranty.^2); 
end
energCin.start_image_number = start_image_number;
energCin.end_image_number = end_image_number;

filename = ['Exp' exp_name '_energCin' ];
outfile = [save_path filename];
save(outfile, 'energCin')

%% Analyse
for i=1:size(energCin.ec,1)
    smoothnEc(i,:)=smoothn(energCin.ec(i,:),10);
end
figure, surf(smoothnEc(:,1:250))


% % Video
% for image_pair_number=start_image_number:end_image_number
%     plot(energCin.ec(:,image_pair_number-start_image_number+1))
%     G(image_pair_number)=getframe;
% end
% movie(G,3,1)