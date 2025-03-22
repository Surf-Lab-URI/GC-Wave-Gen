function L = Lapl(cartDiff,PIVRes,transfo,CST)

% Compute Laplacian based on first order derivatives and extrapolate values
% at the interface. 
% Inputs:
% - cartDiff: structure which contains u_x, u_z, w_x, w_z and Mask
% - PIVRes: structure whic contains x and z PIV coordinates
% - transfo: structure which contains SU (Zeta lines in z coordinates),
%            x_Pix and zeta_pix
% - CST: structure which contains deltaX (m/pix) and deltaT (s/pair)
% 
% Output:
% - L: structure which contains nondimensional Laplacian in x-z,
%      nondimensional flipped Laplacian in x-Zeta and dimensional flipped
%      Laplacian at the surface.

dirX = [0;1]; % directional derivative along x
dirZ = [1;0]; % directional derivative along z
L.u_xx = csapsDiff1(cartDiff.u_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of u ONLY in x direction
L.u_zz = csapsDiff1(cartDiff.u_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of u ONLY in z direction
L.w_xx = csapsDiff1(cartDiff.w_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of w ONLY in x direction
L.w_zz = csapsDiff1(cartDiff.w_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of w ONLY in z direction
L.u_zz = -L.u_zz; % change sign because z is flipped
L.w_zz = -L.w_zz; % change sign because z is flipped
L.Lapl_u = L.u_xx + L.u_zz;
L.Lapl_w = L.w_xx + L.w_zz;
L.Lapl_uInt = TransformVelField_decay_0( flipud(L.Lapl_u.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) );
L.Lapl_wInt = TransformVelField_decay_0( flipud(L.Lapl_w.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) );

% Use the Laplacian extrapolated exactly on the surface
L.u_diff_sum = L.Lapl_uInt(1,:);
L.w_diff_sum = L.Lapl_wInt(1,:);
L.DX = CST.DX;
L.DT = CST.DT;
L.x_pix = transfo.x_pix;
L.zeta_pix = transfo.zeta_pix;

end
