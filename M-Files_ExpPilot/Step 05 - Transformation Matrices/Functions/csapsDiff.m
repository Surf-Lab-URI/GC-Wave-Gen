function [du_dx, du_dz] = csapsDiff(u, smth_param, X, Z)
% This function takes analytic style derivatve. First, it fits a thin shell
% or thin plate through all the data points with a certain amount of rigid-
% ity or smoothing. Then, it find the analytical function of that plate and
% takes the derivatives of the function.
%
%   X = PIVRes.xPIV(3:end-2);
%   Z = PIVRes.zPIV;
%   smth_param = 0.1; % Smoothing Parameter


F = csaps({Z',X}, u, smth_param); % Returns the ppform of a cubic smoothing
% spline to the given data. It fits a thin shell through all the data point
% with a certain amount of rigidity or smoothing. Then, we can calculate or
% find the analytical function of that plate. In fact it finds the function
% that looks like u.
F_INT = fnval( F,{Z',X} ); % Evaluate the function at the Z and X nodes. It
% is a check for the csaps command, and this should be similar to the u.

% Drivatives with smoothing
DF = fndir(F,eye(2)); % Directional derivative of function; it's a function
Diff_F = fnval(DF,{Z',X}); % Evaluates the derivative function.

du_dz = squeeze(Diff_F(1,:,:));
du_dx = squeeze(Diff_F(2,:,:));

end

% The very first way of doing this is the following, but using the analytic
% style derivative is more preferable.

% dudx = zeros(509,985);
% for i = 1:509
%     du = gradient(Cartesian.u(i,:));
%     dx = gradient(PIVRes.xPIV);
%     dudx(i,:) = du./dx;
% end
