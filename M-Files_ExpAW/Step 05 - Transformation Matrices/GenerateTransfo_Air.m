function [transfo] = GenerateTransfo_Air(PixRes, PIVRes, CST)
% This function calculate the transformation matrix of the curvilinear grid
% and at the same time computes the Jacobian. The calculated Jacobian works
% perfectly for x and Zeta coordinates, but not for Xi and Zeta coordinates
% This function is inspired originally from file "transformVelField18.m".


PIVResReal = CST.DX;
PIVDeltaT = CST.DT;
% 1. Interpolate velocities to the curvilinear grid, surface at the bottom,
% and exponentially linearizing surface away from the surface. Perform ffts
% on the LFV.
gravity = 9.81; % m/s^2
gravity_piv = (gravity/PIVResReal) * (PIVDeltaT^2); % Converts the  gravity 
% to pix/(pair^2). Note that the units of pivResReal is m/pix, and the unit
% of pivDeltaT is s/pair.

xx = 1:length(PixRes.PIV_LFV_Surface); % Take the whole LFV surface to extract the longest components
X = find(ismember(PixRes.XPIV_LFV_Surface(xx), PIVRes.xPIV)); % Finds indices that are mem-
% bers of PIV_LFV_Surface and PIVRes.xPIVV; it'll be longer than PIVsurf by a bit to avoid edge effects in SU.
% PixRes.XPIV_LFV_Surface x-location (with respect to PIVAir) of the LFV surface
% PIVRes.xPIV = 4, 8, etc... are the x position of the PIV estimates

surface_fab = 3088-PixRes.PIV_LFV_Surface(xx)+1; %It's a flipud of the PIV_LFV_Surface!   %Fits perfectly with imagesc(flipud(PIV1_A))
surface_fab_smth = 3088-PixRes.PIV_LFV_Surface_smth(xx)+1; % 3088 is the height in pix of a raw PIV Air image 
f = fft(surface_fab);
fs = fft(surface_fab_smth);

%fs = fft(PixRes.Combo_Surface_eta_smth(xx)); 
%f = fft(PixRes.Combo_Surface_eta(xx)); 

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
h = [0 PIVRes.zPIV-CST.GS(end)/2]; % Height above the water surface, i.e. "zeta" in the
% transformations formula. We consider -2 in order to start at the height 2
% instead of 4, which 2 is the mean distance from the surface of a computed
% vector. The 0 is also added as a check.
% OK because X is in pixels res
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
    SU(z_index,:) = h(z_index)+ifft(F.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % zeta+sum(a exp(ikx) exp(-kz))           
    %real(SU) contains location of constant zeta lines to interpolate the velocity to a surface following coordinate system
    SU_E(z_index,:) = h(z_index)+ifft(Fs.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % zeta+sum(a exp(ikx-i om t) exp(-kz))
    %smooth version of SU (to rotate the velocity)
    
    ORB_vel(z_index,:) = ifft(-1i*Fs.*om.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % orbital velocity sum(-i om  a exp(ikx-i om t) exp(-kz))
    ORB_acc(z_index,:) = ifft(-1*Fs.*om.^2.*exp(-h(z_index)*k_LFV),'nonsymmetric'); % orbital acceleration sum(-a om^2 exp(ikx) exp(-kz))
    % smooth version of the orbital velocity
    Jac(z_index,:) =ifft(F.*k_LFV.*exp(-h(z_index)*k_LFV),'nonsymmetric');
    Jac_E(z_index,:) =ifft(Fs.*k_LFV.*exp(-h(z_index)*k_LFV),'nonsymmetric');
end


transfo.SU = SU(:,X); %contains the location in z of constant zeta lines (X is the location of PIV)
transfo.SU_E = SU_E(:,X);
transfo.J = 1 - real(Jac(:,X)); 
transfo.J_smth = 1 - real(Jac_E(:,X)); 
transfo.AK_smth = gradient(real(SU_E(:,X)),PIVRes.GS);
transfo.ORB_vel = ORB_vel(:,X);
transfo.ORB_acc = ORB_acc(:,X);
transfo.zeta_pix=h;
transfo.x_pix=PIVRes.xPIV;


% This Jacobian of transformation is good
% if we just use the x and Zeta coordinates not the Xi and Zeta coordinates
% We have analytically calculated the Jacobian of transformation here which
% means J = 1 - akcos(kx)exp(-kZeta).

% transfo.ORBX = imag(ORB_vel);   % don't erase
% transfo.ORBZ = real(ORB_vel);  % don't erase 
% transfo.ORBX_acc = -imag(ORB_acc);   % don't erase
% transfo.ORBZ_acc = real(ORB_acc);  % don't erase



%% Now SUH which is needed for the inverse horizontal transformation and the scale factors

transfo.SUH = TransformSUField_decay_hor(transfo.SU, PIVRes, imag(transfo.SU));

%% Scale factors

        SUR = real(transfo.SU(2:end, :)); % Real part of the transformation
        % matrix, which contains z while rows are constant Zeta and columns
        % are constant x. That means z(x,zeta).
        SUI = imag(transfo.SU(2:end,:)); % Imaginary part of transformation
        % matrix, which contains xi(x,zeta) - x.
        
        SUHR = real(transfo.SUH(2:end, :)); % Real  part  of transformation
        % matrix, which contains z while rows are constant Zeta and columns
        % are constant Xi. That means z(xi,zeta)  
        
        % Calculating dZ(Xi,Zeta)/dXi and dZ(Xi,Zeta)/dZeta
        Z = smoothn(SUHR, 0, 'robust'); % z = z(xi,zeta)
        Z(~isnan(SUHR)) = SUHR(~isnan(SUHR));
        
        [~, dZ_dZeta] = csapsDiff(Z, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        %dZ_dXi(isnan(SUHR)) = SUHR(isnan(SUHR));
        dZ_dZeta(isnan(SUHR)) = SUHR(isnan(SUHR));
        
        % Calculating dXi(Xi,Zeta)/dx and dXi(Xi,Zeta)/dz 
            X = repmat(PIVRes.xPIV,size(SUR,1),1);
            Xi = SUI + X; % This matrix contains Xi(x,Zeta) , rows are constant
        % Zeta and columns are constant x.
        
            Xi_XZ = inverseTransformVelField_decay(Xi, PIVRes, real(transfo.SU(2:end,:))); % Xi(x,z)
       
        Xi_XZ_Smth = smoothn(Xi_XZ, 0, 'robust'); % Xi = Xi(x,z)
        Xi_XZ_Smth(~isnan(Xi_XZ)) = Xi_XZ(~isnan(Xi_XZ));
        
        [dXi_dx, ~] = csapsDiff(Xi_XZ_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV); % dXi(x,z)/dx, dXi(x,z)/dz
        
            dXi_dx(isnan(Xi_XZ)) = Xi_XZ(isnan(Xi_XZ));
        %dXi_dz(isnan(Xi_XZ)) = Xi_XZ(isnan(Xi_XZ));
        
        % Transform dXi(x,z)/dx, dXi(x,z)/dz to dXi(Xi,Zeta)/dx, dXi(Xi,Zeta)/dz 
            dXi_dX = TransformVelField_decay(dXi_dx, PIVRes, real(transfo.SU(2:end,:)), 'spline');
            dXi_dX = TransformVelField_decay_hor(dXi_dX, PIVRes, imag(transfo.SU(2:end,:)));
        %dXi_dZ = TransformVelField_decay(dXi_dz, PIVRes, real(transfo.SU(2:end,:)));
        %dXi_dZ = TransformVelField_decay_hor(dXi_dZ, PIVRes, imag(transfo.SU(2:end,:)));
        
        % Calculating dX(Xi,Zeta)/dXi and dX(Xi,Zeta)/dZeta 
            X_SF = TransformVelField_decay(X, PIVRes, real(transfo.SU(2:end,:)), 'spline'); % X(x,Zeta)
            X_SFH = TransformVelField_decay_hor(X_SF, PIVRes, imag(transfo.SU(2:end,:))); % X(Xi,Zeta)
        
        X_SFH_Smth = smoothn(X_SFH, 0, 'robust');
        X_SFH_Smth(~isnan(X_SFH)) = X_SFH(~isnan(X_SFH));
        
        [~, dX_dZeta] = csapsDiff(X_SFH_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        %dX_dXi(isnan(X_SFH)) = X_SFH(isnan(X_SFH)); % not as good as 1./dXi_dX
            dX_dZeta(isnan(X_SFH)) = X_SFH(isnan(X_SFH)); % dx/dzeta = - dz/dxi
        
        % To double check, we have to notice that dX_dXi should be equal to
        % the inverse of dXi_dX. 
        
        %% Test 
        % Following calculations is another check for (Xi,Zeta) coordinates
        % Since (Xi,Zeta) are orthogonal then dXi/dXi = 1 and dXi/dZeta = 0
        
        % Xi_SFH = TransformVelField_decay_hor(Xi, PIVRes, imag(transfo.SU(2:end,:)));
        % This matrix contains Xi(Xi,Zeta) while rows are constant Zeta and
        % the columns of matrix are constant Xi.
        
        % Xi_SFH_Smth = smoothn(Xi_SFH, 0, 'robust'); % Xi = Xi(Xi,Zeta)
        % Xi_SFH_Smth(~isnan(Xi_SFH)) = Xi_SFH(~isnan(Xi_SFH));
        
        % [dXi_dXi, dXi_dZeta] = csapsDiff(Xi_SFH_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        % dXi_dXi(isnan(Xi_SFH)) = Xi_SFH(isnan(Xi_SFH));
        % dXi_dZeta(isnan(Xi_SFH)) = Xi_SFH(isnan(Xi_SFH));
        
        %% Scale Factors
        
        dX_dXi = 1./dXi_dX;
        dZ_dXi = -dX_dZeta; % Cauchy-Riemann condition dx/dzeta = - dz/dxi,
        % but dX_dZeta is interpolating better at the edges.
        
        transfo.h1 = sqrt( (dX_dXi).^2 + (dZ_dXi).^2 );
        transfo.h3 = sqrt( (dX_dZeta).^2 + (dZ_dZeta).^2 );
        
        %transfo.Jacobian = (dX_dXi .* dZ_dZeta) - (dX_dZeta .* dZ_dXi); % Jacobian(Xi,Zeta)
        %Jacobian(Xi,Zeta)=h1*h3;


end