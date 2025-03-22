%% Function to compute pressure (in kg/pix/(pair^2))

function [Lapl, press] = compute_pressure(CST,PIVRes,transfo,cartDiff,GrdSpc,Fluid)

%-----------INPUT DATA-----------%
%
% CST: constants and physical parameters
% PIVRes: surface and coordinates in PIV resolution
% transfo: surface acceleration in PIV resolution
% cartDiff: Cartesian derivatives (first order)
% Grdspc: last level of grid space (usually 4 or 16)
% Fluid: if analysis in the air or in water
%
%--------------------------------%

%-----------OUTPUT DATA-----------%
%
% Lapl: Laplacian of Cartesian velocity in PIV resolution
% press: pressure in PIV resolution
%
%---------------------------------%

% Physical parameters
switch Fluid
    case 'WATER'
        % Load file and variables ----------------------------------
        Surf = PIVRes.PF_Surface*GrdSpc;
        x = PIVRes.xPIV;
        U_x = flipud(cartDiff.u_x);
        U_z = flipud(cartDiff.u_z);
        W_z = flipud(cartDiff.w_z);
        W_x = flipud(cartDiff.w_x);
        dx = x(4)-x(3);
        DENSITY = CST.WATER_DENSITY;
        DVISCOSITY = CST.WATER_DVISCOSITY;
        DX = CST.DX_W;
        DT = CST.DT_W;
        z = PIVRes.zPIV;
    case 'AIR'
        % Load file and variables ----------------------------------
        Surf = PIVRes.PF_Surface*GrdSpc;
        x = PIVRes.xPIV;
        U_x = flipud(cartDiff.u_x);
        U_z = flipud(cartDiff.u_z);
        W_z = flipud(cartDiff.w_z);
        W_x = flipud(cartDiff.w_x);
        dx = x(4)-x(3);
        DENSITY = CST.AIR_DENSITY;
        DVISCOSITY = CST.AIR_DVISCOSITY;
        DX = CST.DX;
        DT = CST.DT;
        z = PIVRes.zPIV;
end

% Flags
curved = true;
periodic = false;
Ptop = 0;

% Convert to "pix", "pair"
rhoA_piv = DENSITY*DX^3; % Convert density to kg/(pix^3)
mu_piv = DVISCOSITY*DX*DT; % Convert dynamic viscosity to kg/(pix^3)
Tol = CST.TOLERANCE/DX;

% Compute f from gradients
f = rhoA_piv*(-2*U_x.*W_z + 2*U_z.*W_x);
f = transpose(f);

% Define pressure bottom boundary condition ----------------
% Find surface acceleration

% Smooth acceleration
lambda_cut = 0.01/DX;   % cutoff wavelength in [m] (same units as xLong)
fs = 1/dx;               % sampling wavenumber [cyc/m]
nyq = fs/2;              % Nyquist wavenumber [cyc/m]
wN = 1/(nyq*lambda_cut); % cutoff wavenumber as fraction of Nyquist
[a,b] = butter(2,wN);

% Extract horizontal and vertical acceleration
x_acc = PIVRes.xPIV;
surf_acc = filtfilt(a,b,transfo.ORB_acc(1,:));
utSurf = -imag(surf_acc)';
wtSurf = real(surf_acc)';

% Check that PIV field is correctly identified
if min(abs(x - x_acc) - dx > Tol)
    disp('---> Error: Acceleration surface incorrectly identified')
    % else
    %     disp('Surface acceleration: check')
end

% Viscous term: --------------------------------------------
%%%%%%%%% Insert Laplacian of u (we can calculate directly here)
% Load viscous term data:
dirX = [0;1]; % directional derivative along x
dirZ = [1;0]; % directional derivative along z
Lapl.u_xx = csapsDiff1(cartDiff.u_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of u ONLY in x direction
Lapl.u_zz = csapsDiff1(cartDiff.u_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of u ONLY in z direction
Lapl.w_xx = csapsDiff1(cartDiff.w_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of w ONLY in x direction
Lapl.w_zz = csapsDiff1(cartDiff.w_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of w ONLY in z direction
Lapl.u_zz = -Lapl.u_zz; % change sign because z is flipped
Lapl.w_zz = -Lapl.w_zz; % change sign because z is flipped
Lapl.Lapl_u = Lapl.u_xx + Lapl.u_zz;
Lapl.Lapl_w = Lapl.w_xx + Lapl.w_zz;
switch Fluid
    case 'AIR'
        Lapl.Lapl_uInt = TransformVelField_decay( flipud(Lapl.Lapl_u.*cartDiff.Mask), PIVRes, real(transfo.SU(1:end,:)), 'spline' );
        Lapl.Lapl_wInt = TransformVelField_decay( flipud(Lapl.Lapl_w.*cartDiff.Mask), PIVRes, real(transfo.SU(1:end,:)), 'spline' );
    case 'WATER'
        Lapl.Lapl_uInt = TransformVelField_decay_WATER( Lapl.Lapl_u.*cartDiff.Mask, PIVRes, real(transfo.SU(1:end,:)), 'spline' );
        Lapl.Lapl_wInt = TransformVelField_decay_WATER( Lapl.Lapl_w.*cartDiff.Mask, PIVRes, real(transfo.SU(1:end,:)), 'spline' );
end

% Use the Laplacian extrapolated exactly on the surface
Lapl.u_diff_sum = Lapl.Lapl_uInt(1,:);
Lapl.w_diff_sum = Lapl.Lapl_wInt(1,:);
Lapl.DX = DX;
Lapl.DT = DT;
Lapl.x_pix = transfo.x_pix;
Lapl.zeta_pix = transfo.zeta_pix;

u_diff_sum = squeeze(Lapl.u_diff_sum');
w_diff_sum = squeeze(Lapl.w_diff_sum');

%  Pressure gradient components needed for bottom BC --------
P_x = @(X) interp1(x_acc,-rhoA_piv*utSurf + mu_piv*u_diff_sum,X,'spline','extrap');
P_z = @(Z) interp1(x_acc,-rhoA_piv*wtSurf + mu_piv*w_diff_sum,Z,'spline','extrap');

% Smoothed version of the pressure gradient
P_x = @(X) interp1(x_acc,smooth(-rhoA_piv*utSurf + mu_piv*u_diff_sum,0.05),X,'spline','extrap');
P_z = @(Z) interp1(x_acc,smooth(-rhoA_piv*wtSurf + mu_piv*w_diff_sum,0.05),Z,'spline','extrap');

switch Fluid
    case 'WATER'
        % Create pressure field ------------------------------------
        Out = solve_Poisson_WATER(f,Surf,x,z,P_x,P_z,curved,periodic,Ptop);
        % unpack the output:
        p2 = Out.P;       % pressure field [ kg/pix/(pair^2) ]
        % fields: 'P','z','x','surf','surf_index'

        % Return p to original dimensions
        p = nan(length(x),length(z));
        pz_num = size(p2,2);
        p(:,length(z) - pz_num+1:length(z)) = p2;

        % Save as mat file -------------------------
        press.p = flipud(p');
        press.x = x;
        press.z = z;
    case 'AIR'
        % Create pressure field ------------------------------------
        Out = solve_Poisson(f,Surf,x,z,P_x,P_z,curved,periodic,Ptop);
        % unpack the output:
        p2 = Out.P;       % pressure field [ kg/pix/(pair^2) ]
        % fields: 'P','z','x','surf','surf_index'

        % Return p to original dimensions
        p = nan(length(x),length(z));
        pz_num = size(p2,2);
        p(:,length(z) - pz_num+1:length(z)) = p2;

        % Save as mat file -------------------------
        press.p = flipud(p');
        press.x = x;
        press.z = z;
end
