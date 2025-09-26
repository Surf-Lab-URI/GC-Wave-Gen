tic,
clear all
close all
%% Paramètres d'entrée
exp_name = 'LC2_dt7ms_1';
piv_res = 40d-6;
deltaT = 7d-3; %s
GrdSpc = [128 64 32 16 8];
Res=piv_res;
start_image_number = 250;
end_image_number = 450;

%% Paramètres intrinsèques
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];
%% 
for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    %%
    u=compVel.delx_ints*Res/deltaT;
    v=compVel.dely_ints*Res/deltaT;
    % Calcul des champs de vitesses dans le repere de surface
    u_rs = interpTransfo255472(compVel.surf1,u,compVel.SU,GrdSpc);
    v_rs = interpTransfo255472(compVel.surf1,v,compVel.SU,GrdSpc);
    % Moyenne selon x dans le repere de surface
    um_rs=nanmean(u_rs,2);
    vm_rs=nanmean(v_rs,2);
    % Creation d'une matrice a la taille du probleme
    um_rs_mat=repmat(um_rs,1,size(u,2));
    vm_rs_mat=repmat(vm_rs,1,size(v,2));
    % Matrice turbulente dans le repere de surface
    up_rs=u_rs-um_rs_mat;
    vp_rs=v_rs-vm_rs_mat;
    % Matrice turbulente dans le repere d'origine
    up = interpTransfo255472Inverse(compVel.surf1,up_rs,compVel.SU,GrdSpc);
    vp = interpTransfo255472Inverse(compVel.surf1,vp_rs,compVel.SU,GrdSpc);
    % Calcul des differents
%     up2m_rs=nanmean(up_rs.^2,2);
%     vp2m_rs=nanmean(vp_rs.^2,2);
%     upvpm_rs=nanmean(up_rs.*vp_rs,2);
%     ecpm_rs=nanmean((up_rs.^2+vp_rs.^2)/2,2);

    ReyDec.up=up;
    ReyDec.vp=vp;
    ReyDec.up_rs=up_rs;
    ReyDec.vp_rs=vp_rs;
    ReyDec.um_rs=um_rs;
    ReyDec.vm_rs=vm_rs;
%     ReyDec.up2m_rs=up2m_rs;
%     ReyDec.vp2m_rs=vp2m_rs;
%     ReyDec.upvpm_rs=upvpm_rs;
%     ReyDec.ecpm_rs=ecpm_rs;
    
    filename = ['Exp' exp_name '_ReyDec_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'ReyDec')
  toc  

end
