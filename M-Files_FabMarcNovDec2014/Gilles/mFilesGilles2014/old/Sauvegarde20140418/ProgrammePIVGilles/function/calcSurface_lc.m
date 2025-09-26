function Surface = calcSurface_lc(surfbig,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_right)

s2 = surfbig.z_s;
imgsurf=surfbig.img;
s2 = size(imgsurf,1) - s2;
% figure, plot(s2)
s2=fliplr(s2);

GS=GrdSpc(end);
IW=IntrWndw(end);
%% function compute orbitals
%%  calculate orbitals
%VELOCITY RECONSTRUCTION
% xx = 0:1:13100-1;
% s2 = sin(pi*xx/5000);
dx = pivRes;
x_s = lim_left_big+GS:GS:lim_left_big+lim_right-1-GS;
z_s = 0:GS:2048-GS*2;
SU = zeros(length(z_s),length(x_s));
f = fft(s2); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s2)-1)]/(length(s2)-1)*2*pi/dx; %wavenumber of mode

for j = 1:length(z_s)
    su=0;
    for i = 1:floor(length(k)/2)%/10
        su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i));
    end
        
    SU(j,:)=su-z_s(j);
  end

Surface.SU=(SU);
% g=0;
% for i=1:floor(length(s2)/2)
% g=g+(fa(i))*cos(k(i)*x_s(1:end-1)+fp(i));
% end
% g=g-fa(1)/2;
% g=interp1(x_s(1:end-1),g,x_s,'linear','extrap');%X_s(2:end)

