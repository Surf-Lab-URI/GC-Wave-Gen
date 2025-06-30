%% Spatial Grid
clear
clc
close all

U_s = 0;
U_inf = 5;
h_a = 5e-4;
h_w = 1e-2;

dz = 1e-5;
z = (0:dz:h_a*40)';

iS = find(z == 0);

nZ = length(z);

D2 = zeros(nZ);

for i = 2:length(z)-1
    D2(i, (i-1):(i+1)) = [1,-2,1]/dz^2;
end

U = zeros(nZ,1);

U(1:iS) = 0;%U_s*exp(z(1:iS)/h_w);
U(iS:end) = 0;%U_inf - (U_inf - U_s)*exp(-z(iS:end)/h_a);


plot(U,z)
%%
k = 2*pi/15;

A = diag(U)*D2-k^2*diag(U)-diag(D2*U);
B = -k^2*eye(nZ)+D2;

[V,D] = eig(A,B);
%%
c = diag(D);
real(c);   
imag(c);

figure
plot(abs(V(:,22)),z);
hold on
plot(angle(V(:,22)),z);

%% Spatial Grid V2
U_s = 0;
U_inf = 5;
h_a = 5e-4;

dz = 1e-5;
z = (0:dz:h_a*40)';

iS = find(z == 0);

nZ = length(z);

D2 = zeros(nZ);

for i = 2:length(z)-1
    D2(i, (i-1):(i+1)) = [1,-2,1]/dz^2;
end

U = zeros(nZ,1);

U(1:iS) = U_s;%U_s*exp(z(1:iS)/h_w);
U(iS:end) = 0;%U_inf - (U_inf - U_s)*exp(-z(iS:end)/h_a);


A = diag(U)*D2-k^2*diag(U)-diag(D2*U);
B = -k^2*eye(nZ)+D2;

[V,D] = eig(A,B);



%% Chebyshev Spectral Method with discretized integration: This doesn't work unless you use super super huge nX
clear
clc
close all

N = 20;

nX = 100000;
x = linspace(-0.999999999,0.999999999,nX)';
dx = x(2)-x(1);
U = 1;

c = ones(N+1,1);
c(1) = 2;

b = zeros(N+1,1);

parfor n = 0:N
    b(n+1) = 2./(pi*c(n+1)).*dx.*sum(U.*chebyshevT(n,x).*(1-x.^2).^(-0.5));
end

U_r = zeros(nX,1);
for n = 0:N
    U_r = U_r + b(n+1)*chebyshevT(n,x);
end

figure
plot(U_r(20:80))
%% With symbolic integration
clear
clc
close all

N = 20;

c = ones(N+1,1);
c(1) = 2;

syms x

U = x^2+1;

b = zeros(N+1,1);

parfor n = 0:N
    b(n+1) = vpa(2/pi/c(n+1)*int(U*chebyshevT(n,x)*(1-x.^2).^(-0.5),x,-1,1))
end

U_r = 0;
for n = 0:N
    U_r = U_r + b(n+1)*chebyshevT(n,x);
end
U_r

%% Calculate coefficient matrix contributions from terms with U multiplied by eigenfunction or its derivatives

m = 1;

A = zeros(1,N+1);

for p = 0:N
    for q = 0:(N-1)
        if p + q == m
            for r = (q+1):N
                if mod(r + q,2) == 1
                    A(1,r) = A(1,r) + 1/2*b(p+1)*2/c(q+1)*r;
                end
            end
        end
    end
end