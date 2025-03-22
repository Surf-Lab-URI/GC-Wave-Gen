function [uInt] = TransformVelField_decay_0(u, PIVRes, SU)
% This function is vertically transforming z to the zeta, e.g. the velocity
% u(x,z) will be transformed to the u(x,zeta). It is inspired from: MFiles\
% transformation\transformVelField18.m.
%
%   Notice that the transformation matrix should be considered from  second
%   row, since the first row is zeta = 0, i.e. the surface exactly.
%   u = flipud(Cartesian.u.*Cartesian.Mask);
%   w = flipud(Cartesian.u.*Cartesian.Mask);
%   SU = real(transfo.SU(2:end,:));
% 
% %-% Added velocity extrapolated at z = 0. First line corresponds to surface velocity 
% %-% (Fabio 01/28/2022)

GS = PIVRes.GS;
SU = (length(PIVRes.zPIV)-(7799-SU+1)/GS) +1; %length(PIVRes.zPIV) == size(u,1); % correction by Fabrice
% % SU=(length(PIVRes.zPIV)-(7799-SU+1)/GS) +1; %length(PIVRes.zPIV) == size(u,1);
% % % (7799-SU+1) is right side up
% % % (7799-SU+1)/GS is right side up in PIV 
% % % length(PIVRes.zPIV)-(7799-SU+1)/GS) +1 is up side down in PIV

%%% This was used for compute_pressure, gives the same results at the surface:
%SU = (SU-1)/GS; % Change SU from pixel to PIV resolution. We Substract -1
% because of flipping; the same as what we did for "PIVRes.PF_Surface = (pi
% xRes.PF_Surface(PIVRes.xPIV) - 1)/CompVel.GS" in compute velocity.
%%%

PFSurf = PIVRes.PF_Surface ; 
%PFSurf = size(u,1)-PIVRes.PIVFused_Surface+1; %flip upside down

D = [];
uInt = NaN(size(u,1)+1,size(u,2));

for col = 1:size(SU,2)
    
    ggu = u(:,col);
    fvp = find(~isnan(ggu), 1, 'first'); % First Velocity Position: Returns
    % the first element that is not NaN, which is the surface.
    
    if isempty(fvp)
        continue
    end
    
    ggu = ggu(fvp:end); % The velocity from the surface to the end.
    
    zi = fvp - PFSurf(col) : size(u,1) - PFSurf(col); % Initial Data Sites: 
    % The location of velocity measurments in the z-direction. 
    hh = SU(:,col);
    zt = hh(hh < size(u,1) + 50) - PFSurf(col); % Target data sites. Notice
    % that the SU was NEVER flipped. We added 50 to be sure that we capture
    % everything.
        
%     % Start interpolating at z = 0
    zt = [0;zt];

    %D = [D (fvp - PFSurf(col))];% Checks whether the location of first PIV
    % measurement is between 0 and 1. We can also check this 
    %D = [D zt(1)]; That's the average target height measurement (should average to 0.5)
    
    % We only extrapolate on the bottom (near the surface) , not on the top
    % (away from surface).
    ggu_Int1 = interp1(zi,ggu,zt, 'spline', 'extrap');
    ggu_Int2 = interp1(zi,ggu,zt, 'linear', NaN); % For the linear methods,
    % interp1 return NaN for any element of zt that is outside the interval
    % spanned by zi. For other methods or when it is specified by 'extrap',
    % interp1 performs extrapolation for out of range values.
    ggu_Int = ggu_Int1;
    ggu_Int(find(~isnan(ggu_Int2), 1, 'last') + 1:end) = NaN;
    uInt(1:length(ggu_Int),col) = ggu_Int;
    
end

% if min(D) < 0
%     display 'There are some measurements below the surface!!!';
% end

end
