function Transfo = calcTransfo(surfbig,pivRes,lim_left_big,lim_right)
%% Remise de la surface aux dimensions du probleme
s2 = surfbig.z_s;
s2 = size(surfbig.img,1) - s2;
s2=fliplr(s2);
s_pivsurf=s2;
%% Parametres
dx = pivRes;
z_s = 0:1:2048-1;
x_s = lim_left_big+1:1:lim_left_big+lim_right;
%% Decomposition de Fourier
f = fft(s_pivsurf); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s_pivsurf)-1)]/(length(s_pivsurf)-1)*2*pi/dx; %wavenumber of mode
%% Reconstruction de la surface avec une decroissance exponentielle liée au nombre d'onde k
for j = 1:length(z_s)
    g = 0;
    h = 0;
    su=0;
    for i = 1:1:floor(length(k)/2/1)% Choix du nomdre de "k" que l'on soufaite garder
        su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); 
    end
    SU(j,:)=su-z_s(j); 
end
Transfo=(SU);



