

%Z=compVel.xPIV;
%X=compVel.zPIV;
%u= compVel.INTdelx;
function [dudx, dudz] = csapsDiff(u, smth_param, X, Z)
% smth_param = 0.1; littlesmoothing
% smth_param = 0.001; some smoothing
% smth_param = 1d-5; lots smoothing
% smth_param = 1-1d-4; no smoothing
% analytic style derivatve
%
F = csaps({Z',X}, u, smth_param); % finds the function that looks like u
%
F_INT=fnval(F,{Z',X}); % evaluate the function at the {} nodes: interpolation
%
%derivatives with smoothing
DF = fndir(F,eye(2)); % derivative of the function; it's a function too
Diff_F = fnval(DF,{Z',X}); % evaluated the derivative function
dudz = squeeze(Diff_F(1,:,:));
dudx = squeeze(Diff_F(2,:,:));



%[cartDiff.u_x, cartDiff.u_z] = csapsDiff(cart.u(:,3:end-2), 0.001, pivRes.xPIV(3:end-2), pivRes.zPIV);
%    [cartDiff.w_x, cartDiff.w_z] = csapsDiff(cart.w(:,3:end-2), 0.001, pivRes.xPIV(3:end-2), pivRes.zPIV);
%    cartDiff.u_z = -cartDiff.u_z;
%    cartDiff.w_z = -cartDiff.w_z;
