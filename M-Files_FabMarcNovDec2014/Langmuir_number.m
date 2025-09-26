% calculate Langmuir number
function La=Langmuir_number(a,k,U10)

    
    T=0.073; % T is the surface tension coefficient
    kv=1e-6; % kv is the
%     u_star=1.1e-3; %u_star is the friction velocity (1994 Komen)
    u_star_w=sqrt(1.2/1000*1.1d-3*U10^2);
    g=9.8;       %gravity accelaration
    ak=sqrt(2)*a*k;
%     DX=5.65571e-05;
%     ak=sqrt(2)*std(a*DX)*k;
    rho=1000;

    sigma=sqrt(g.*k.*(1+T.*k.^2/rho/g));
    down=(sigma./k).*(ak).^2.*u_star_w^2;
    La=(kv.*k).^3/down;
%     Laa=kv^3*k^2/(sigma*a*a*u_star^2) if here a^2 is a*a, then the same.
    La=sqrt(La);
 end


