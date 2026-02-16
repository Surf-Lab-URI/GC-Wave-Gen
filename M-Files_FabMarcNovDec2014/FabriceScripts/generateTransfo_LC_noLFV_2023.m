function transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes)
%
% inspired from E:\MFiles\transformation\transformVelField18.m
pivDeltaT = compVel.DT;%10d-6;  
pivResReal = compVel.DX;%5.65d-7; %m/pix
% load('E:\data\Exp4\ComputedVelocities\Exp4_compVel_1004.mat')
%%
altitude=[0 pivRes.zPIV-pivRes.GS/2]; %compVel.GS/2 on average with moving surface
X=pivRes.GS:pivRes.GS:2048-pivRes.GS;
%% 1. interpolate velocities to curvilinear grid, surface at the bottom, and exponentially linearizing surface away from surface
% do ffts on lfv
gravity = 9.81; %m/s^2
gravity_piv = gravity/pivResReal*pivDeltaT^2;
%
% X = find(ismember(pixRes.xLFV, pivRes.xPIV))-1;
%X = 1:length(Surface_PIV);

%
f  = fft(Surface_PIV); %Fourier Modal decomposition
%% build sym_vec, vector that takes care of symmetry and scales amplitude of fft
sym_vec = nan(size(f));
sym_vec(1) = 1; sym_vec(length(f)/2+1) = 1; 
sym_vec(2:length(f)/2) = 2;
sym_vec(length(f)/2+2:end) = 0;
%
F = f.*sym_vec;
N = length(f);
%%
k_lfv = (0:(N-1))/(N)*2*pi; %wavenumber of mode
%
% h     = [0 2:4:altitude]; % 
h     = altitude;
%OK because X is in pixels res
%
%%
SU_E = nan(length(h),length(X));
% AK_E = nan(length(h),length(X));
ORB_E = nan(length(h),length(X));
% ORB_E_acc = nan(length(h),length(X));
%
%loop on heights above water surface
for z_index = 1:length(h)
    su_e = 0; orb_e = 0;% orb_e_acc = 0;
    % loop on wavenumbers
    for i = 1:(N/2+1)/4 % limit to modes all the way to PIV resolution GS, no need to go beyong 
        expD = exp(-h(z_index)*k_lfv(i)); % exponential decay
        Ealpha = 1i .* k_lfv(i);
        E    =  exp( Ealpha .* (X-1) ); % cos(kx)
%         dE = Ealpha .* E;
        ORB_SCALE = sqrt( gravity_piv * k_lfv(i) ); % omega
        %
        su_e =  su_e + expD .* F(i) .* E; 
%         ak_e =  ak_e + expD .* F(i) .* dE; 
        orb_e = orb_e + expD .* F(i) .* E .* ORB_SCALE;
%         orb_e_acc = orb_e_acc + expD .* F(i) .* E .* ORB_SCALE.^2;
    end
    %
    SU_E(z_index,:) = su_e./N + h(z_index);
    %
%     AK_E(z_index,:) = ak_e./N; % only 10 first modes
    ORB_E(z_index,:) = orb_e./N; % only 10 first modes
%     ORB_E_acc(z_index,:) = -orb_e_acc./N; % only 10 first modes
end
%
transfo.SU = real(SU_E);
%
% transfo.AK = real(AK_E);
% transfo.ORB = ORB_E;
% transfo.ORB_acc = ORB_E_acc;
%
transfo.ORBX = -real(ORB_E);  % don't erase
transfo.ORBZ = -imag(ORB_E);   % don't erase
% transfo.ORBX_acc = -real(ORB_acc);  % don't erase
% transfo.ORBZ_acc = imag(ORB_acc);   % don't erase

%
end