clear all
close all
dx = 1/4000;
L = 1-1/4000;
x = 0:dx:L;
left = 1001;
right = 3000;
%% lfv high res
k = [5 2 2 50];
a = [1 1 5 0.3];
for i=1:length(k)
    z(i,:) = a(i)*sin(2*pi*k(i)*x);
end
z_lfv_hr = sum(z,1);
figure, plot(z_lfv_hr)
%
%% pivsurf high res
x_ps = x(left:right);
z_ps_hr = z_lfv_hr(left:right);
figure, plot(z_ps_hr)
%% lfv low res
k_lr = [5 2 2];
a_lr = [1 1 5];
for i=1:length(k_lr)
    z_lr(i,:) = a_lr(i)*sin(2*pi*k_lr(i)*x);
end
z_lfv_lr = sum(z_lr,1);
figure, plot(z_lfv_lr)
%
%% do ffts on lfv low res
f = fft(z_lfv_lr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(z_lfv_lr)-1)]/(length(z_lfv_lr)-1)*2*pi/dx; %wavenumber of mode
h = 0:1:1000;

for j = 1:length(h)
%     g = 0;
%     h = 0;
    su=0;
    for i = 1:50%floor(length(k)/2)
        su=su+fa(i)*exp(-h(j)*k(i))*cos(k(i)*x(1:end)+fp(i)); % first 5 modes for surface
    end
    SU(j,:)=su-h(j);
end

figure, plot(z_lfv_lr,'k', 'LineWidth', 5), hold on, plot(SU(1,:),'g' , 'LineWidth', 2)
clear SU

%% do ffts on pivsurf high res
f = fft(z_ps_hr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(z_ps_hr)-1)]/(length(z_ps_hr)-1)*2*pi/dx; %wavenumber of mode
% h = 0:1:1000;
for j = 1:length(h)
    su=0;
    for i = 1:50%floor(length(k)/2)
        su=su+fa(i)*exp(-h(j)*k(i))*cos(k(i)*x(1:length(z_ps_hr))+fp(i)); % first 5 modes for surface
    end
    SU_ps(j,:)=su-h(j);
end
figure, plot(z_ps_hr,'k', 'LineWidth', 5), hold on, plot(SU_ps(1,:),'g' , 'LineWidth', 2)
clear SU

%% do ffts on lfv hr
f = fft(z_lfv_hr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(z_lfv_hr)-1)]/(length(z_lfv_hr)-1)*2*pi/dx; %wavenumber of mode
h = 0:1:1000;

for j = 1:length(h)
%     g = 0;
%     h = 0;
    su=0;
    for i = 1:100%floor(length(k)/2)
        su=su+fa(i)*exp(-h(j)*k(i))*cos(k(i)*x(1:end)+fp(i)); % first 5 modes for surface
    end
    SU(j,:)=su-h(j);
end

figure, plot(z_lfv_hr,'k', 'LineWidth', 5), hold on, plot(SU(1,:),'g' , 'LineWidth', 2)
clear SU

%% patchwork
z_patch = [z_lfv_lr(1:left-1) z_lfv_hr(left:right) z_lfv_lr(right+1:end)];
figure, plot(z_patch)

%% do ffts on z_patch
f = fft(z_patch); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k = [0:(length(z_patch)-1)]/(length(z_patch)-1)*2*pi/dx; %wavenumber of mode
h = 0:1:1000;

for j = 1:length(h)
%     g = 0;
%     h = 0;
    su=0;
    for i = 1:100%floor(length(k)/2)
        su=su+fa(i)*exp(-h(j)*k(i))*cos(k(i)*x(1:end)+fp(i)); % first 5 modes for surface
    end
    SU(j,:)=su-h(j);
end

figure, plot(z_patch,'k', 'LineWidth', 5), hold on, plot(SU(1,:),'g' , 'LineWidth', 2)
clear SU

%% rebuild pivsurf high res from low modes from lfv + high modes from pivsurf

%% do ffts on lfv
f = fft(z_lfv_lr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k_lfv = [0:(length(z_lfv_lr)-1)]/(length(z_lfv_lr)-1)*2*pi/dx; %wavenumber of mode
h = 0:1:1000;

for j = 1:length(h)
%     g = 0;
%     h = 0;
    su=0;
    for i = 1:3%floor(length(k_lfv)/2)
        su=su+fa(i)*exp(-h(j)*k_lfv(i))*cos(k_lfv(i)*x(left:right)+fp(i)); % first 5 modes for surface
    end
    SU(j,:) = su;
end

% SU = SU(:,left:right);
% hold on, plot(SU(1,:),'r')

f = fft(z_ps_hr); %Fourier Modal decomposition
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f); %phase of mode
k_ps = [0:(length(z_ps_hr)-1)]/(length(z_ps_hr)-1)*2*pi/dx; %wavenumber of mode

for j = 1:length(h)
%     g = 0;
%     h = 0;
    su=0;
    for i = 3:50%floor(length(k_ps)/2)
        su=su+fa(i)*exp(-h(j)*k_ps(i))*cos(k_ps(i)*x(1:length(z_ps_hr))+fp(i)); % first 5 modes for surface
    end
    SU(j,:) = SU(j,:) + su - h(j);
end

figure, plot(z_ps_hr,'k', 'LineWidth', 5), hold on, plot(SU(1,:),'g' , 'LineWidth', 2)









