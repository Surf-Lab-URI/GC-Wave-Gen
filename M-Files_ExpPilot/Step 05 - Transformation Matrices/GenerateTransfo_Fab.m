function [transfo] = GenerateTransfo_Fab(PixRes, PIVRes, CST)
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

xx=[6539:15572]; %500 points on each side of the PIV images
%xx=[1:length(PixRes.Combo_Surface_X)];
X = find(ismember(PixRes.Combo_Surface_X(xx), PIVRes.xPIV)); % Finds indices that are mem-
% bers of Combo_Surface_X and PIVRes.xPIVV; it'll be longer than PIVsurf by a bit to avoid edge effects in SU.
% PixRes.Combo_Surface_X x location (with respect to PIVFused) of the combo surf
% PIVRes.xPIV = 4, 8, etc... are the x position of teh PIV estimates

surface_kianoosh=3034-PixRes.Combo_Surface_eta-2; %surface I need to do SU with to be
%like Kianoosh; It's a flipud of the Combo_Surface!
surface_kianoosh_smth=3034-PixRes.Combo_Surface_eta_smth-1; %3034 is the height in pix of a raw PIV image 
f = fft(surface_kianoosh(xx));
fs = fft(surface_kianoosh_smth(xx));

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
h = [0 PIVRes.zPIV-2]; % Height above the water surface, i.e. "zeta" in the
% transformations formula. We consider -2 in order to start at the height 2
% instead of 4, which 2 is the mean distance from the surface of a computed
% vector. The 0 is also added as a check.
% OK because X is in pixels res
%initialize 
SU = nan(length(h),length(X));
SU_E = nan(length(h),length(X));
ORB_E = nan(length(h),length(X));
Jac = nan(length(h),length(X));
Jac_E = nan(length(h),length(X));
for z_index = 1:length(h) % Loop on heights above the water surface
    
    su = 0; jac = 0;
    su_e = 0; orb_e = 0; jac_e = 0;
%     orb_e_acc = 0; 
    
    for i = 1:N/(2*PIVRes.GS/2)+1 % Loop on the wave numbers dont do all the modes; only down to PIV resolution
        expD = exp(-h(z_index)*k_LFV(i)); % Exponentially decaying term
        Ealpha = 1i .* k_LFV(i);
        E = exp( Ealpha .* (X-1) ); % cos(kx)
        ORB_SCALE = sqrt( gravity_piv * k_LFV(i) ); % Omega
        
        su = su + expD .* F(i) .* E;
        su_e =  su_e + expD .* Fs(i) .* E;
%         orb_e = orb_e + expD .* Fs(i) .* E .* ORB_SCALE;
        orb_e = orb_e + expD .* F(i) .* E .* ORB_SCALE;
        jac = jac + expD .* F(i) .* k_LFV(i) .* E;
        jac_e = jac_e + expD .* Fs(i) .* k_LFV(i) .* E;
    end
    
   
    SU(z_index,:) = su./N + h(z_index);
    SU_E(z_index,:) = su_e./N + h(z_index);
    ORB_E(z_index,:) = orb_e./N; % Only 10 first modes
%     ORB_E_Acc(z_index,:) = -orb_e_acc./N; % Only 10 first modes
    Jac(z_index,:) = jac./N;
    Jac_E(z_index,:) = jac_e./N;
end

transfo.SU = SU; %contains the location in z of constant zeta lines (X is the location of PIV)
transfo.SU_E = SU_E;
transfo.J = 1 - real(Jac); 
transfo.J_smth = 1 - real(Jac_E); 
transfo.AK_smth = gradient(real(SU_E),PIVRes.GS);
transfo.ORB_smth = ORB_E;
transfo.zeta_pix=h;
transfo.x_pix=PIVRes.xPIV;

% This Jacobian of transformation is good
% if we just use the x and Zeta coordinates not the Xi and Zeta coordinates
% We have analytically calculated the Jacobian of transformation here which
% means J = 1 - akcos(kx)exp(-kZeta).

% transfo.ORBX = -real(ORB);  % don't erase
% transfo.ORBZ = imag(ORB);   % don't erase
% transfo.ORBX_acc = -real(ORB_acc);  % don't erase
% transfo.ORBZ_acc = imag(ORB_acc);   % don't erase


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
            dXi_dX = TransformVelField_decay(dXi_dx, PIVRes, real(transfo.SU(2:end,:)));
            dXi_dX = TransformVelField_decay_hor(dXi_dX, PIVRes, imag(transfo.SU(2:end,:)));
        %dXi_dZ = TransformVelField_decay(dXi_dz, PIVRes, real(transfo.SU(2:end,:)));
        %dXi_dZ = TransformVelField_decay_hor(dXi_dZ, PIVRes, imag(transfo.SU(2:end,:)));
        
        % Calculating dX(Xi,Zeta)/dXi and dX(Xi,Zeta)/dZeta 
            X_SF = TransformVelField_decay(X, PIVRes, real(transfo.SU(2:end,:))); % X(x,Zeta)
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