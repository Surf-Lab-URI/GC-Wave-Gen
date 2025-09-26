function SU = SUcc(surfLfv, surff, ResPiv, hResLfv, vResLfv, GS, Att)

%% % % Lfv % % %%
%% Remettre la surface a l'echelle
s = surfLfv.z_s;
s1 = s*vResLfv/ResPiv; % surface in PIVpix
xi = 0:hResLfv:(length(s)-1)*hResLfv;  % initial data sites
xt = 0:ResPiv:(length(s)-1)*hResLfv; % target data sites
% xt = 0:(length(s)-1)*hResLfv; % target data sites

ss2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)
s2 = size(surfLfv.img,1)*hResLfv - ss2;
s2 = s2 - mean(s2);
s_lfv=s2;
clear s s1 xi xt ss2 s2
%% Extraire le coefficient de decroissance dominant
f = fft(s_lfv); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s_lfv)-1)]/(length(s_lfv)-1)*2*pi; %wavenumber of mode
[C,I] = max(fa);
kdom=k(I);
clear f fa fp k C
%% % % Pivsurfcc % % %%
s = surff.z_s_f; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surfcc = size(surff.img,1) - s;   
% x_s=GS:GS:size(surf1.img,2)-GS;
x_s=GS:GS:size(surff.img,2)-GS;
z_s=0:GS:size(surff.img,1)-2*GS;
if GS==1 % On veut que l'image d'entree et de sortie soit a la meme taille 
    x_s=1:GS:size(surff.img,2);
    z_s=0:GS:size(surff.img,1)-GS;
end

f = fft(surfcc); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fa(1)=fa(1)/2;
fp = angle(f); %phase of mode
k = [0:(length(surfcc)-1)]/(length(surfcc)-1)*2*pi;
for j = 1:length(z_s)
    su=0;
    for i = 1:floor(length(k)/2/GS)
        su=su+fa(i)*exp(-Att*z_s(j)*kdom)*cos(k(i)*x_s(1:end)+fp(i));
    end
    SU(j,:)=su-z_s(j);
end

