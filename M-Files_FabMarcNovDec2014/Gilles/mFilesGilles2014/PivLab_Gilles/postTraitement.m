
%% Parametres
% Filtre par ecart-type
stdthresh = 7;
% Filtre
epsilon = 0.1;
thresh = 5;

load('E:\ComputedVelocities\20140804\LC1_1\ExpLC1_1_compVel_PivLab_0368.mat')
%%
u = utable; v= vtable;
u=inpaint_nans(u,4); v=inpaint_nans(v,4); % On interpole les NaN
% Filtre median local
[u,v] = medfiltCheck(u,v,epsilon,thresh);
u=inpaint_nans(u,4); v=inpaint_nans(v,4);
% Filtre globale par ecart-type
u = stdDevCheck(u,stdthresh); v = stdDevCheck(v,stdthresh);
u=inpaint_nans(u,4); v=inpaint_nans(v,4);
% Leger smootn du resultat
u=smoothn(u,0.5); v=smoothn(v,0.5);
