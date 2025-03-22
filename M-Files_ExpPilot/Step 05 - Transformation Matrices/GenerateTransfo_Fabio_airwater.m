function [Uorb,Vorb] = GenerateTransfo_Fabio_airwater(XPIVW_PIVSurfW1_Surface, PIVW_PIVSurfW1_Surface, CST,h2, w, PIV1_W)
% This function calculate the transformation matrix of the curvilinear grid
% and at the same time computes the Jacobian. The calculated Jacobian works
% perfectly for x and Zeta coordinates, but not for Xi and Zeta coordinates
% This function is inspired originally from file "transformVelField18.m".


PIVResReal = CST.DX_W;
PIVDeltaT = CST.DT_W;
PIVRes.xPIV=CST.IW/2:CST.GS:(w-CST.IW/2);
PIVRes.zPIV=CST.IW/2:CST.GS:(h2-CST.IW/2);
PIVRes.GS = CST.GS;
if mod(length(PIVW_PIVSurfW1_Surface),2) == 1
    PIVW_PIVSurfW1_Surface = PIVW_PIVSurfW1_Surface(2:end);
    XPIVW_PIVSurfW1_Surface = XPIVW_PIVSurfW1_Surface(2:end);
end

b=ones(64,1)/64;
PIVW_PIVSurfW1_Surface_smth=filtfilt(b,1,PIVW_PIVSurfW1_Surface);

% 1. Interpolate velocities to the curvilinear grid, surface at the bottom,
% and exponentially linearizing surface away from the surface. Perform ffts
% on the LFV.
gravity = 9.81; % m/s^2
gravity_piv = (gravity/PIVResReal) * (PIVDeltaT^2); % Converts the  gravity 
% to pix/(pair^2). Note that the units of pivResReal is m/pix, and the unit
% of pivDeltaT is s/pair.

I1 = find(XPIVW_PIVSurfW1_Surface==1);
Iend = find(XPIVW_PIVSurfW1_Surface==size(PIV1_W,2));
X = I1:Iend; % We save only the surface of PIVSurf Water that overlaps with PIV Water
% xx = find(ismember(X, PIVRes.xPIV+I1-1));
xx = 1:length(X);

surface_fab=3091-PIVW_PIVSurfW1_Surface+1; %It's a flipud of the Combo_Surface!   %Fits perfectly with imagesc(flipud(FusedPIV1))
surface_fab_smth=3091-PIVW_PIVSurfW1_Surface_smth+1; %3091 is the height in pix of a raw PIV image 
f = fft(surface_fab);
fs = fft(surface_fab_smth);

% Build Symmetric Vector: The vector that takes care of symmetry and scales
% amplitude of the fft. Sym_Vec = [1 2 2 ... 2 2 1 0 0 ... 0]
Sym_Vec = NaN(size(f));
Sym_Vec(1) = 1; Sym_Vec(length(f)/2+1) = 1; 
Sym_Vec(2:length(f)/2) = 2;
Sym_Vec(length(f)/2+2:end) = 0;

F = f.*Sym_Vec; 
Fs = fs.*Sym_Vec;
N = length(f);


k_LFV = ( (0:(N-1))/N ) * (2*pi); % Wavenumber of mode or Fourier transform
% h = [0 PIVRes.zPIV-CST.GS(end)/2]; % Height above the water surface, i.e. "zeta" in the
% transformations formula. We consider -2 in order to start at the height 2
% instead of 4, which 2 is the mean distance from the surface of a computed
% vector. The 0 is also added as a check.
% OK because X is in pixels res
h = 0:h2;
eta = surface_fab;
om = nan(1,length(eta)); % Omega
om(1:N/2+1) = sqrt(gravity_piv * k_LFV(1:N/2+1));
om(N/2+2:N) = om(N/2:-1:2);

%inititalize 
SU = nan(length(h),N);
SU_E = nan(length(h),N);
ORB_vel = nan(length(h),N);
ORB_acc = nan(length(h),N);
Jac = nan(length(h),N);
Jac_E = nan(length(h),N);

for z_index = 1:length(h)
    SU(z_index,:) = -h(z_index)+ifft(F.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % zeta+sum(a exp(ikx) exp(-kz))           
    %real(SU) contains location of constant zeta lines to interpolate the velocity to a surface following coordinate system
    SU_E(z_index,:) = -h(z_index)+ifft(Fs.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % zeta+sum(a exp(ikx-i om t) exp(-kz))
    %smooth version of SU (to rotate the velocity)
    
    ORB_vel(z_index,:) = ifft(-1i*Fs.*om.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % orbital velocity sum(-i om  a exp(ikx-i om t) exp(-kz))
    ORB_acc(z_index,:) = ifft(-1*Fs.*om.^2.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % orbital acceleration sum(-a om^2 exp(ikx) exp(-kz))
    % smooth version of the orbital velocity
    Jac(z_index,:) =ifft(F.*k_LFV.*exp(-h(z_index)*k_LFV),'nonsymmetric');
    Jac_E(z_index,:) =ifft(Fs.*k_LFV.*exp(-h(z_index)*k_LFV),'nonsymmetric');
end


transfo.SU = SU(:,X(xx)); %contains the location in z of constant zeta lines (X is the location of PIV)
transfo.SU_E = SU_E(:,X(xx));
transfo.J = 1 - real(Jac(:,X(xx))); 
transfo.J_smth = 1 - real(Jac_E(:,X(xx))); 
transfo.AK_smth = gradient(real(SU_E(:,X(xx))),PIVRes.GS);
transfo.ORB_vel = ORB_vel(:,X(xx));
transfo.ORB_acc = ORB_acc(:,X(xx));
transfo.zeta_pix = h;

% This Jacobian of transformation is good
% if we just use the x and Zeta coordinates not the Xi and Zeta coordinates
% We have analytically calculated the Jacobian of transformation here which
% means J = 1 - akcos(kx)exp(-kZeta).

% transfo.ORBX = imag(ORB_vel);   % don't erase
% transfo.ORBZ = real(ORB_vel);  % don't erase 
% transfo.ORBX_acc = -imag(ORB_acc);   % don't erase
% transfo.ORBZ_acc = real(ORB_acc);  % don't erase

ORB_vel_water = (transfo.ORB_vel(2:end,:));
PF_Surface = PIVW_PIVSurfW1_Surface(X);
PIVRes.PF_Surface = size(PIV1_W,1)-PF_Surface+1;
PIVRes.GS = 1;
PIVRes.zPIV = 1:size(ORB_vel_water,1);
PIVRes.xPIV = 1:size(ORB_vel_water,2);
Uorb = inverseTransformVelField_decay_WATER(-imag(ORB_vel_water), PIVRes, real(transfo.SU(2:end,:))); 
Vorb = inverseTransformVelField_decay_WATER(real(ORB_vel_water), PIVRes, real(transfo.SU(2:end,:))); 

end