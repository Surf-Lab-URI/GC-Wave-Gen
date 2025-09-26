%% Evolution de l'énergie cinetique totale et de l'energie cinetique due aux vitesses orbitales en fonction de la profondeur (definie a partir de SU) et du temps. 
tic,
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC1_dt7ms_1';
% Res = 1.2700e-04; %m/pixel A VERIFIER (depend de la vue CC ou longitudinale)
piv_res = 40d-6;
deltaT = 7d-3; %s
GrdSpc = [128 64 32 16 8];
Res=piv_res*GrdSpc(end);
start_image_number = 300;
end_image_number = 700;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% Calcul de l'energie cinetique 
for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
%% Energie cinetique totale au temps image_pair_number
    % Remet la surface a zero en utilisant la transformation SU
    delx_tranfo = interpTransfo255472(compVel.surf1,compVel.delx_ints,compVel.SU,GrdSpc);
    dely_tranfo = interpTransfo255472(compVel.surf1,compVel.dely_ints,compVel.SU,GrdSpc);
    % Calcul la moyenne en fonction de la profondeur des composantes
    % longitudianale et verticale de la vitesse totale.
    delx_mean=nanmean(delx_tranfo,2);
    dely_mean=nanmean(dely_tranfo,2);
%% Energie cinetique orbitale au temps image_pair_number
    delxOrb_tranfo = interpTransfo255472(compVel.surf1,compVel.delxOrb,compVel.SU,GrdSpc);
    delyOrb_tranfo = interpTransfo255472(compVel.surf1,compVel.delyOrb,compVel.SU,GrdSpc);
    % Calcul la moyenne en fonction de la profondeur des composantes
    % longitudianale et verticale des vitesses orbitales
    delxOrb_mean=nanmean(delxOrb_tranfo,2);
    delyOrb_mean=nanmean(delyOrb_tranfo,2);
    
    delxTurb = delx_tranfo - delxOrb_tranfo;
    delyTurb = dely_tranfo - delyOrb_tranfo;
    delxTurb_mean = nanmean(delxTurb,2);
    delyTurb_mean = nanmean(delyTurb,2);

    %% Sauvegarde des informations calculee au temps image_pair_number
    EcAndEcOrb.delx_mean(:,image_pair_number-start_image_number+1)=delx_mean;
    EcAndEcOrb.dely_mean(:,image_pair_number-start_image_number+1)=dely_mean;
    EcAndEcOrb.delxOrb_mean(:,image_pair_number-start_image_number+1)=delxOrb_mean;
    EcAndEcOrb.delyOrb_mean(:,image_pair_number-start_image_number+1)=delyOrb_mean;
    EcAndEcOrb.delxTurb_mean(:,image_pair_number-start_image_number+1)=delxTurb_mean;
    EcAndEcOrb.delyTurb_mean(:,image_pair_number-start_image_number+1)=delyTurb_mean;
    EcAndEcOrb.ec(:,image_pair_number-start_image_number+1)=1/2*(delx_mean.^2+dely_mean.^2);
    EcAndEcOrb.ecOrb(:,image_pair_number-start_image_number+1)=1/2*(delxOrb_mean.^2+delyOrb_mean.^2);
end

EcAndEcOrb.start_image_number = start_image_number;
EcAndEcOrb.end_image_number = end_image_number;

%% Sauvegarde
filename = ['Exp' exp_name '_calcEcOrb2' ];
outfile = [save_path filename];
save(outfile, 'EcAndEcOrb')
toc
%% Analyse
% for i=1:size( EcAndEcOrb.ec,1)
%     smoothnEc(i,:)=smoothn(EcAndEcOrb.ec(i,:),10);
%     smoothnEcOrb(i,:)=smoothn(EcAndEcOrb.ecOrb(i,:),10);
% end
% figure, surf(smoothnEc-smoothnEcOrb)





