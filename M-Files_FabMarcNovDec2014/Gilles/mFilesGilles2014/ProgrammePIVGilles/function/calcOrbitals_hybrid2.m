function Orbitals = calcOrbitals_hybrid2(surfLfv,surfbig,vRes,hRes,pivRes,deltaT,pivResReal,IntrWndw,GrdSpc,lim_left_big,lim_left,lim_right)

s = surfLfv.z_s;
s1 = s*vRes; % surface in PIVpix
xi = 0:hRes:(length(s)-1)*hRes;  % initial data sites
yi = 0:vRes:(size(surfLfv.img,1)-1)*vRes;
xt = 0:pivRes:(length(s)-1)*hRes; % target data sites
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



%%

    s = surfLfv.z_s;
    s = size(surfLfv.img,1) - s;
    s = s - mean(s);
    b=ones(1,50)/50;
    s=filtfilt(b,1,s);
    s2 = surfbig.z_s;
    s2 = size(surfbig.img,1) - s2;
    s2=fliplr(s2);
    s2 = s2 - mean(s2);
    
    a=s(1:4434);
    b=s(4435:4435+5286-1); bm=s2;
    c=s(4434+5286:end); 
    a=a-mean(b(1:150))+mean(bm(1:150));
    b_left=b-mean(b(1:150))+mean(bm(1:150));
    c=c-mean(b(end-150:end))+mean(bm(end-150:end));   
    b_right=b-mean(b(end-150:end))+mean(bm(end-150:end));
    tz = 50;
    bm(1:tz) = 1/(tz-1)*((0:tz-1).*bm(1:tz) + (tz-1:-1:0).*b_left(1:tz));
    bm(end-tz+1:end) = 1/(tz-1)*((0:tz-1).*b_right(end-tz+1:end) + (tz-1:-1:0).*bm(end-tz+1:end));
    t=[a bm c];
    
    lim_left=5127;
    dx = pivRes;
    x_s = lim_left+GS:GS:lim_left+lim_right-1-GS;
    z_s = 0:GS:2048-GS*2;
    u = zeros(length(z_s),length(x_s));
    w = zeros(length(z_s),length(x_s));
    f = fft(t); %Fourier Modal decomposition
    fa = 2*abs(f)/length(f); %amplitude of mode
    fp = angle(f); %phase of mode
    k = [0:(length(t)-1)]/(length(t)-1)*2*pi/dx; %wavenumber of mode

    for j = 1:length(z_s)
        g = 0;
        h = 0;
        su=0;
           for i = 1:floor(length(k)/2/GS)
            g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only low modes for orbital vel from lfv_surf
            h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
            su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i));
           end
         u(j,:) = g;  % combine low and hig modes for the orbital vel.
         w(j,:) = h;
         SU(j,:)=su-z_s(j);
    end
    

Orbitals.u=u;
Orbitals.w=w;
Orbitals.SU=(SU);

