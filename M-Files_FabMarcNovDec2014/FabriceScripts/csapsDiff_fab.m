

%Z=compVel.xPIV;
%X=compVel.zPIV;
%u= compVel.INTdelx;
function [dudx, dudz] = csapsDiff_fab(u, smth_param, DX)
% smth_param = 0.1; littlesmoothing
% smth_param = 0.001; some smoothing
% smth_param = 1d-5; lots smoothing
% smth_param = 1-1d-4; no smoothing
% analytic style derivatve
%
% DX is the actual DX in the same units as u (for example, it's m if u is
% in m/s

Z=[1:size(u,1)]; % THIS assume dx=dz=1 for the calculation below;
X=[1:size(u,2)];

F = csaps({Z',X}, u, smth_param); % finds the function that looks like u
%
F_INT=fnval(F,{Z',X}); % evaluate the function at the {} nodes: interpolation
%
%derivatives with smoothing
DF = fndir(F,eye(2)); % derivative of the function; it's a function too
Diff_F = fnval(DF,{Z',X}); % evaluated the derivative function
dudz = squeeze(Diff_F(1,:,:));
dudx = squeeze(Diff_F(2,:,:));

%ATTENTION!!!! YOU NEED to divide dudx by the "true" dx after this opration
DZ=DX; % assume same in both directions
dudz=dudz/DZ;
dudx=dudx/DX;


%[cartDiff.u_x, cartDiff.u_z] = csapsDiff(cart.u(:,3:end-2), 0.001, pivRes.xPIV(3:end-2), pivRes.zPIV);
%    [cartDiff.w_x, cartDiff.w_z] = csapsDiff(cart.w(:,3:end-2), 0.001, pivRes.xPIV(3:end-2), pivRes.zPIV);
%    cartDiff.u_z = -cartDiff.u_z;
%    cartDiff.w_z = -cartDiff.w_z;
