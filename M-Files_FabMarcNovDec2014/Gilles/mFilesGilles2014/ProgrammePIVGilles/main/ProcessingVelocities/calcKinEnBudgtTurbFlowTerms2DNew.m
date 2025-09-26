tic,
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC2_dt7ms_1';
piv_res = 40d-6;
deltaT = 7d-3; %s
deltaPiv = 1/7.2;
GrdSpc = [128 64 32 16 8];
Res=piv_res*GrdSpc(end);
start_image_number = 251;
end_image_number = 449;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
nu=1e-6; % Voir les donnees des thermistors pour une estimation plus precise


%% N.B. Pertinence de ces calculs pour une turbulence #D non isotropique (Langmuir) en ne possedant qu'un champ de vitesse 2D?


for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number-1)]);
    U(:,1)=ReyDec.um_rs; % U(i,j,t), t=1==n-1, t=2 ==n, t=3==n+1
    V(:,1)=ReyDec.vm_rs;
    u(:,:,1)=ReyDec.up_rs;
    v(:,:,1)=ReyDec.vp_rs;
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
    U(:,2)=ReyDec.um_rs;
    V(:,2)=ReyDec.vm_rs;
    u(:,:,2)=ReyDec.up_rs;
    v(:,:,2)=ReyDec.vp_rs;
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number+1)]);
    U(:,3)=ReyDec.um_rs;
    V(:,3)=ReyDec.vm_rs;
    u(:,:,3)=ReyDec.up_rs;
    v(:,:,3)=ReyDec.vp_rs;
    clear ReyDec
    %% Evolution (schema centre en temps)
    Evol = (1/2)*(nanmean(u(:,:,3).^2,2)-nanmean(u(:,:,1).^2,2)+ nanmean(v(:,:,3).^2,2)-nanmean(v(:,:,1).^2,2))/(2*deltaPiv);
    %% Advection (schema centre en espace) - entre le terme dans la direction manquante et le terme terme moyenne en espace, il ne reste plus qu'un terme...

    for j = 2:size(v,2)-1
        dpx(:,j-1)=(u(2:end-1,j+1,2).^2-u(2:end-1,j-1,2).^2+v(2:end-1,j+1,2).^2-v(2:end-1,j-1,2).^2)/(2*piv_res*GrdSpc(end));
    end
    for i = 2:size(v,1)-1
        dpy(i-1,:)=(u(i+1,2:end-1,2).^2-u(i-1,2:end-1,2).^2+v(i+1,2:end-1,2).^2-v(i-1,2:end-1,2).^2)/(2*piv_res*GrdSpc(end));
    end
    Advec = (1/2)*(U(2:end-1,2).*nanmean(dpx,2)+V(2:end-1,2).*nanmean(dpy,2));
    clear dpx dpy
    
    %% Turbulent transport of (schema centre en espace) - idem
    % pression - inconnu, peut etre calculer a partir de l'equqtion de l'energie cinetique turbulente et des autres termes. Voir le residu si apres. 
    % energie cinetique
    for j = 2:size(v,2)-1
        dpx(:,j-1)=( (u(2:end-1,j+1,2).^2+v(2:end-1,j+1,2).^2).*u(2:end-1,j+1,2) - (u(2:end-1,j-1,2).^2+v(2:end-1,j-1,2).^2).*u(2:end-1,j-1,2) ) / (2*piv_res*GrdSpc(end));
    end
    
    for i = 2:size(v,1)-1
        dpy(i-1,:)=( (u(i+1,2:end-1,2).^2+v(i+1,2:end-1,2).^2).*v(i+1,2:end-1,2) - (u(i-1,2:end-1,2).^2+v(i-1,2:end-1,2).^2).*v(i-1,2:end-1,2) ) / (2*piv_res*GrdSpc(end));
    end
    TurbTransport = -(1/2)*(nanmean(dpx,2)+nanmean(dpy,2));
    clear dpx dpy

    % Molecular viscous transport
    for j = 2:size(v,2)-1
        dudx(:,j-1) = (u(2:end-1,j+1,2)-u(2:end-1,j-1,2)) / (2*piv_res*GrdSpc(end));
        dvdx(:,j-1) = (v(2:end-1,j+1,2)-v(2:end-1,j-1,2)) / (2*piv_res*GrdSpc(end));
    end
    for i = 2:size(v,1)-1
        dudy(i-1,:) = (u(i+1,2:end-1,2)-u(i-1,2:end-1,2)) / (2*piv_res*GrdSpc(end));
        dvdy(i-1,:) = (v(i+1,2:end-1,2)-v(i-1,2:end-1,2)) / (2*piv_res*GrdSpc(end));
    end

    for j = 3:size(v,2)-2 % Les indices pour "u" et "dudx" sont decales de 1
        dpx(:,j-2) = (u(3:end-2,j+1,2).*dudx(2:end-1,j) + (1/2)*v(3:end-2,j+1,2).*(dudy(2:end-1,j)+dvdx(2:end-1,j)) - u(3:end-2,j-1,2).*dudx(2:end-1,j-2) - (1/2)*v(3:end-2,j-1,2).*(dudy(2:end-1,j-2)+dvdx(2:end-1,j-2))) / (2*piv_res*GrdSpc(end));
    end

    for i = 3:size(v,1)-2 % Les indices pour "u" et "dudx" sont decales de 1
        dpy(i-2,:) = (v(i+1,3:end-2,2).*dvdy(i,2:end-1) + (1/2)*u(i+1,3:end-2,2).*(dudy(i,2:end-1)+dvdx(i,2:end-1)) - u(i-1,3:end-2,2).*dudx(i-2,2:end-1) - (1/2)*v(i-1,3:end-2,2).*(dudy(i-2,2:end-1)+dvdx(i-2,2:end-1))) / (2*piv_res*GrdSpc(end));
    end
    
    MolVisquTransport = 2*nu*(nanmean(dpx,2)+nanmean(dpy,2));
    
    clear dpx dpy
    %% Production
    for i=2:size(v,1)-1
        Prod(i-1) = -(nanmean(v(i,:,2).^2,2).*((V(i+1,2)-V(i-1,2))/(2*piv_res*GrdSpc(end))) + nanmean(u(i,:,2).*v(i,:,2),2).*((U(i+1,2)-U(i-1,2))/(2*piv_res*GrdSpc(end))) );
    end
    Prodd = Prod';
    %% Dissipation -
    Dissip3D=-6*nu*nanmean(dudx.^2+dudy.^2+dudy.*dvdx,2);
    Dissip2D=-nu*nanmean(2*dudx.^2+2*dvdy.^2+(dvdx+dudy).^2,2);
    

    %% Residu
    Residu=Evol(3:end-2)+Advec(2:end-1)-TurbTransport(2:end-1)-MolVisquTransport-Prodd(2:end-1)-Dissip3D(2:end-1);
    
%     KEBTF.Evol=Evol;
%     KEBTF.Advec=Advec;
%     KEBTF.TurbTransport=TurbTransport;
%     KEBTF.Prod=Prodd;
%     KEBTF.MolVisquTransport = MolVisquTransport;
%     KEBTF.Dissip=Dissip3D;
%     KEBTF.Residu=Residu;
    
    KEB.Evol(:,image_pair_number-start_image_number+1) = Evol(3:end-2);
    KEB.Advec(:,image_pair_number-start_image_number+1) = Advec(2:end-1);
    KEB.TurbTransport(:,image_pair_number-start_image_number+1) = TurbTransport(2:end-1);
    KEB.Prod(:,image_pair_number-start_image_number+1) = Prodd(2:end-1);
    KEB.MolVisquTransport(:,image_pair_number-start_image_number+1) = MolVisquTransport;
    KEB.Dissip(:,image_pair_number-start_image_number+1) = Dissip3D(2:end-1);
    KEB.Residu(:,image_pair_number-start_image_number+1) = Residu;

    
    
%     % Sauvegarde
%     filename = ['Exp' exp_name '_KinEnBudgtTurbFlow_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
%     outfile = [save_path filename];
%     save(outfile, 'KEBTF')
end
toc

KEB.start_image_number = start_image_number
KEB.end_image_number = end_image_number

% Sauvegarde
filename = ['Exp' exp_name '_KinEnBudgtTurbFlow'];
outfile = [save_path filename];
save(outfile, 'KEB')

% figure, plot(Evol(3:end-2))
% hold on, plot(Advec(2:end-1),'r')
% hold on, plot(TurbTransport(2:end-1),'g')
% hold on, plot(MolVisquTransport,'m')
% hold on, plot(Prod(2:end-1),'k')
% hold on, plot(Dissip3D(2:end-1),'c')
% hold on, plot(Residu,'o')

