clear all
close all
dx = 1/4000;
L = 1-1/4000;
x = 0:dx:L;
left = 1001;
right = 3000;
%
%% lfv high res
k = [5 2 2 50];
a = [1 1 5 0.3];
for i=1:length(k)
    z(i,:) = a(i)*sin(2*pi*k(i)*x);
end
z_lfv_hr = sum(z,1);
% figure, plot(z_lfv_hr)
%
%% pivsurf high res
x_ps = x(left:right);
z_ps_hr = z_lfv_hr(left:right);
% figure, plot(z_ps_hr)
%
%% do ffts on pivsurf
f = fft(z_ps_hr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(z_ps_hr)-1)]/(length(z_ps_hr)-1)*2*pi/dx; %wavenumber of mode
h = 0:1:1000;
for j = 1:length(h)
    su=0;
    for i = 1:50%floor(length(k)/2)
        su_temp = fa(i)*exp(-h(j)*k(i))*cos(k(i)*x(1:length(z_ps_hr))+fp(i));
        su = su + su_temp;
%         keyboard;
    end
    SU_ps(j,:)=su-h(j);
end
figure, plot(z_ps_hr), hold on, plot(SU_ps(1,1:2000),'r')