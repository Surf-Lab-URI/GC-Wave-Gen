function Orbitals = calcOrbitals_hybrid(surfLfv,surfbig,vRes,hRes,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right)

s = surfLfv.z_s;
s1 = s*vRes; % surface in PIVpix
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(surfLfv.img,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
ss2 = spline(xi,s1,xt);  % smth_surface_rescaled_vert_splined (spline interpolation)
s2 = size(surfLfv.img,1)*hRes - ss2;
s2 = s2 - mean(s2);
s_lfv=s2;

s2 = surfbig.z_s;
s2 = size(surfbig.img,1) - s2;
% figure, plot(s2)
s2=fliplr(s2);
s_pivsurf=s2;


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
       for i = 1:11%floor(length(k)/2)%/10
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only low modes for orbital vel from lfv_surf
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
       end
    u(j,:) = g;
    w(j,:) = h;
end


dx = pivRes;
x_s = lim_left_big+GS:GS:lim_left_big+lim_right-1-GS;
f = fft(s_pivsurf); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s_pivsurf)-1)]/(length(s_pivsurf)-1)*2*pi/dx; %wavenumber of mode

for j = 1:length(z_s)
    g = 0;
    h = 0;
    su=0;
    for i = 1:5
        su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % first 5 modes for surface
    end
    for i = 6:floor(length(k)/2/GS)%/10 Cut to the resolution of the final PIV (i.e. 1/GS), no need to go to higher k
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only hig modes for orbital vel from piv_surf
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
        su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % rest of modes for surface
    end
    u(j,:) = u(j,:)+g;  % combine low and hig modes for the orbital vel.
    w(j,:) = w(j,:)+h;
    SU(j,:)=su-z_s(j);
end

Orbitals.u=u;
Orbitals.w=w;
Orbitals.SU=(SU);



