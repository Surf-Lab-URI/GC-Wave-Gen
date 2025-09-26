tic,
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC2_dt7ms_3';
piv_res = 40d-6;
deltaT = 7d-3; %s
GrdSpc = [128 64 32 16 8];
Res=piv_res*GrdSpc(end);
start_image_number = 002;
end_image_number = 299;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
nu=1e-6; % Voir les donnees des thermistors pour une estimation plus precise


%% N.B. Pertinence de ces calculs pour une turbulence #D non isotropique (Langmuir) en ne possedant qu'un champ de vitesse 2D?


for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number-1)]);
    VarM=ReyDec;
    VarM.up2m_rs=nanmean(VarM.up_rs.^2,2);
    VarM.vp2m_rs=nanmean(VarM.vp_rs.^2,2);
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
    Var=ReyDec;
    Var.up2m_rs=nanmean(Var.up_rs.^2,2);
    Var.vp2m_rs=nanmean(Var.vp_rs.^2,2);
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number+1)]);
    VarP=ReyDec;
    VarP.up2m_rs=nanmean(VarP.up_rs.^2,2);
    VarP.vp2m_rs=nanmean(VarP.vp_rs.^2,2);
    clear ReyDec
    %% Evolution (schema centre en temps)
    Evol = (((VarP.up2m_rs+VarP.vp2m_rs)-(VarM.up2m_rs+VarM.vp2m_rs))/(2*2*deltaT))';
    %% Advection (schema centre en espace) - entre le terme dans la direction manquante et le terme terme moyenne en espace, il ne reste plus qu'un terme...
    for i = 2:size(Var.vm_rs)-1
        Advec(i) = Var.vm_rs(i)*((Var.up2m_rs(i+1,:)+Var.vp2m_rs(i+1,:))-(Var.up2m_rs(i-1,:)+Var.vp2m_rs(i-1,:)))/(2*2*piv_res*GrdSpc(end));
    end
    %% Turbulent transport of (schema centre en espace) - idem
    % pression - inconnu, peut etre calculer a partir de l'equqtion de l'energie cinetique turbulente et des autres termes. Voir le residu si apres. 
    % energie cinetique
    TurbTransport=0;
    for i = 2:size(Var.vm_rs)-1
        TurbTransport(i) = (nanmean((Var.up_rs(i+1,:).^2+Var.vp_rs(i+1,:).^2).*Var.vp_rs(i+1,:),2)-nanmean((Var.up_rs(i-1,:).^2+Var.vp_rs(i-1,:).^2).*Var.vp_rs(i-1,:),2))/(2*2*piv_res*GrdSpc(end));
    end
    % 
    
    
    for i = 2:size(Var.vm_rs)-1
        for j = 2:size(Var.vm_rs)-1
    QttAMoy(i,j)=(Var.up_rs(i,j)+Var.vp_rs(i,j)).*((Var.up_rs(i+1,j)-Var.up_rs(i-1,j))+(Var.vp_rs(i,j+1)-Var.vp_rs(i,j-1)))/(2*2*piv_res*GrdSpc(end));
        end
        TurbTransport(i)=TurbTransport(i)+nu*nanmean(QttAMoy(i,:),2);
    end
    
    
    
    %% Production
    for i=2:size(Var.vm_rs)-1
        Prod(i)=nanmean(Var.up_rs(i,:).* Var.vp_rs(i,:).*((Var.um_rs(i+1)-Var.um_rs(i-1))/(2*piv_res*GrdSpc(end))),2);
    end
    %% Dissipation -
    for i=2:size(Var.vm_rs)-1
        Dissip(i)=15*nu*(nanmean(((Var.vm_rs(i+1,:)-Var.vm_rs(i-1,:))/(2*piv_res*GrdSpc(end))).^2,2));
    end
    
    
    
    
    %% Residu
    Residu(2:size(Var.vm_rs)-1)=Evol(2:size(Var.vm_rs)-1)+Advec(2:size(Var.vm_rs)-1)-TurbTransport(2:size(Var.vm_rs)-1)-Prod(2:size(Var.vm_rs)-1)-Dissip(2:size(Var.vm_rs)-1);
    
    KEBTF.Evol=Evol;
    KEBTF.Advec=Advec;
    KEBTF.TurbTransport=TurbTransport;
    KEBTF.Prod=Prod;
    KEBTF.Dissip=Dissip;
    KEBTF.Residu=Residu;
    
    %% Sauvegarde
    filename = ['Exp' exp_name '_KinEnBudgtTurbFlow_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'KEBTF')
end


% figure, plot(Evol)
% hold on, plot(Advec,'r')
% hold on, plot(TurbTransport,'g')
% hold on, plot(Prod,'k')
% hold on, plot(Dissip,'c')
% hold on, plot(Residu,'m')

