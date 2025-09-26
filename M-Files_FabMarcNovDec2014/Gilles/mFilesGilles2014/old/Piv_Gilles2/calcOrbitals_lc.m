function Orbitals = calcOrbitals_lc(surf,pivRes,iws,deltaT,pivResReal)

s2 = surf;
gravity=9.81; %m/s^2
gravity_piv=gravity/pivResReal*deltaT^2;

%% function compute orbitals
%%  calculate orbitals
%VELOCITY RECONSTRUCTION
% xx = 0:1:13100-1;
% s2 = sin(pi*xx/5000);
dx = pivRes;
x_s = pivRes * iws : pivRes * iws : (length(s2)-1)*pivRes;
z_s = pivRes * iws : pivRes * iws : 2048-pivRes * iws; %vertical coordinate up to 20cm with 1mm resolution (to be adjusted to match PIV)
u = zeros(length(z_s),length(x_s)-1);
w = zeros(length(z_s),length(x_s)-1);
f = fft(s2); %Fourier Modal decompositionfigure
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx; %wavenumber of mode
for j = 1:length(z_s)
    g = 0;
    h = 0;
    for i = 1:floor(length(s2)/2)%/10
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i));
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i));
    end
    u(j,:) = g;
    w(j,:) = h;
end
Orbitals.u=u;
Orbitals.w=w;
% g=0;
% for i=1:floor(length(s2)/2)
% g=g+(fa(i))*cos(k(i)*x_s(1:end-1)+fp(i));
% end
% g=g-fa(1)/2;
% g=interp1(x_s(1:end-1),g,x_s,'linear','extrap');%X_s(2:end)

