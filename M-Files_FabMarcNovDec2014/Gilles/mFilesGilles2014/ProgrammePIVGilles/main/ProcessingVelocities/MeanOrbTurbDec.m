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
    % Vitesse totale
    u=compVel.delx_ints*Res/deltaT;
    v=compVel.dely_ints*Res/deltaT;
    % Vitesse orbitale
    uorb=compVel.delxOrb*Res/deltaT;
    vorb=compVel.delyOrb*Res/deltaT;
   
    u_rs=interpTransfo255472(compVel.surf1,u,compVel.SU,GrdSpc);
    v_rs=interpTransfo255472(compVel.surf1,v,compVel.SU,GrdSpc);
    
    uorb_rs=interpTransfo255472(compVel.surf1,uorb,compVel.SU,GrdSpc);
    vorb_rs=interpTransfo255472(compVel.surf1,vorb,compVel.SU,GrdSpc);
    
    uMoinsUorb_rs=u_rs-uorb_rs;
    vMoinsVorb_rs=v_rs-vorb_rs;
    
    uMoinsUorbm_rs=nanmean(uMoinsUorb_rs,2);
    vMoinsVorbm_rs=nanmean(vMoinsVorb_rs,2);
    
    uMoinsUorbm_rs_mat=repmat(uMoinsUorbm_rs,1,size(u,2));
    vMoinsVorbm_rs_mat=repmat(vMoinsVorbm_rs,1,size(v,2));
 
    up_rs=u_rs-uorb_rs-uMoinsUorbm_rs_mat;
    vp_rs=v_rs-vorb_rs-vMoinsVorbm_rs_mat;
    
    up = interpTransfo255472Inverse(compVel.surf1,up_rs,compVel.SU,GrdSpc);
    vp = interpTransfo255472Inverse(compVel.surf1,vp_rs,compVel.SU,GrdSpc);
    
    up2m_rs=nanmean(up.^2,2);
    vp2m_rs=nanmean(vp.^2,2);
    upvpm_rs=nanmean(up.*vp,2);
    ecpm_rs=nanmean((up.^2+vp.^2)/2,2);
    
    MOTDec.up=up;
    MOTDec.vp=vp;
    MOTDec.up_rs=up_rs;
    MOTDec.vp_rs=vp_rs;
    MOTDec.uMoinsUorbm_rs = uMoinsUorbm_rs;
    MOTDec.vMoinsVorbm_rs = vMoinsVorbm_rs;
    MOTDec.up2m_rs=up2m_rs;
    MOTDec.vp2m_rs=vp2m_rs;
    MOTDec.upvpm_rs=upvpm_rs;
    MOTDec.ecpm_rs=ecpm_rs;
    
    filename = ['Exp' exp_name '_MOTDec' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
    outfile = [save_path filename];
    save(outfile, 'MOTDec')
   
    
toc
end