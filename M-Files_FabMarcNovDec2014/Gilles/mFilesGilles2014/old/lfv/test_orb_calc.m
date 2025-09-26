clear all
close all

%% compute orbitals on large field of view
dx = 1;
Lv = 200;
L = 1000;
x_s = 0:dx:L-1;
s = 3*sin(pi*x_s/200) + sin(pi*x_s/50);

figure, plot(s)

gravity = 9.81;

z_s = 0:1:Lv; %vertical coordinate up to 20cm with 1mm desolution (to be adjusted to match PIV)
u = zeros(length(z_s),length(x_s)-1);
w = zeros(length(z_s),length(x_s)-1);
f = fft(s); %Fourrier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s)-1)]/(length(s)-1)*2*pi/dx; %wavenumber of mode
for j = 1:length(z_s)
    g = 0;
    h = 0;
%     j
    for i=1:floor(length(s)/2)%/10
        g = g-fa(i)*sqrt(gravity*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i));%+mean(s)/floor(length(s))/2*sqrt(gravity*k(i))*exp(-z_s(j)*k(i));
        h = h+fa(i)*sqrt(gravity*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i));%-mean(s)/floor(length(s))/2*sqrt(gravity*k(i))*exp(-z_s(j)*k(i));
    end
    u(j,:) = g;
    w(j,:) = h;
end
u1 = u;
w1 = w;
figure, imagesc(u1(2:end,:)), colorbar

%% subsample

u11 = u(:,400:600);
figure, imagesc(u11(2:end,:)), colorbar
w11 = w(:,400:600);
s11 = s(400:600);
figure, plot(s11);

%% compute orbitals on subsample of surface
dx = 1;
L = 200;
x_s = 0:dx:L-1;
% s = 3*sin(pi*x_s/200) + sin(pi*x_s/50);
s = s11;
figure, plot(s)

gravity = 9.81;
%%

z_s = 0:1:Lv; %vertical coordinate up to 20cm with 1mm desolution (to be adjusted to match PIV)
u = zeros(length(z_s),length(x_s)-1);
w = zeros(length(z_s),length(x_s)-1);
f = fft(s); %Fourrier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(s)-1)]/(length(s)-1)*2*pi/dx; %wavenumber of mode
for j = 1:length(z_s)
    g = 0;
    h = 0;
%     j
    for i=1:floor(length(s)/2)%/10
        g = g-fa(i)*sqrt(gravity*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i));%+mean(s)/floor(length(s))/2*sqrt(gravity*k(i))*exp(-z_s(j)*k(i));
        h = h+fa(i)*sqrt(gravity*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i));%-mean(s)/floor(length(s))/2*sqrt(gravity*k(i))*exp(-z_s(j)*k(i));
    end
    u(j,:) = g;
    w(j,:) = h;
end
u2 = u;
w2 = w;
figure, imagesc(u2(2:end,:)), colorbar

