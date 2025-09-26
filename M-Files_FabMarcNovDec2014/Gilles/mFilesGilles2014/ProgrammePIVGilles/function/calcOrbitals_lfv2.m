function Orbitals = calcOrbitals_lfv2(surfLfv,surfbig,vRes,hRes,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right)

s = surfLfv.z_s; %From LFV
s1 = s*vRes; % surface in PIVpix
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(surfLfv.img,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
ss2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)
s2 = size(surfLfv.img,1)*hRes - ss2;
s2 = s2 - mean(s2);
s_lfv=s2;




gravity=9.81; %m/s^2
gravity_piv=gravity/pivResReal*deltaT^2;

GS=GrdSpc(end);
IW=IntrWndw(end);


%% function compute orbitals
%%  calculate orbitals
%VELOCITY RECONSTRUCTION
% xx = 0:1:13100-1;
% s2 = sin(pi*xx/5000);
dx = pivRes;
x_s = lim_left+GS:GS:lim_left+lim_right-1-GS;
z_s = 0:GS:2048-GS*2;
u = zeros(length(z_s),length(x_s));
w = zeros(length(z_s),length(x_s));
f = fft(s_lfv); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s_lfv)-1)]/(length(s_lfv)-1)*2*pi/dx; %wavenumber of mode

for j = 1:length(z_s)
    g = 0;
    h = 0;
       for i = 1:floor(length(k)/2/GS)%/10
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only low modes for orbital vel from lfv_surf
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
       end
    u(j,:) = g;
    w(j,:) = h;
end



Orbitals.u=u;
Orbitals.w=w;




