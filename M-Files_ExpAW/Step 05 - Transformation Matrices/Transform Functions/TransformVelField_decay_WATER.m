function [uInt] = TransformVelField_decay_WATER(u, PIVRes, SU, Type)
% This function is vertically transforming z to the zeta, e.g. the velocity
% u(x,z) will be transformed to the u(x,zeta). It is inspired from: MFiles\
% transformation\transformVelField18.m.
%
%   Notice that the transformation matrix should be considered from  second
%   row, since the first row is zeta = 0, i.e. the surface exactly.
%   But if we want uINT at the surface, we should consider the first row.
% 
%   u = flipud(Cartesian.u.*Cartesian.Mask);
%   w = flipud(Cartesian.w.*Cartesian.Mask);
%   SU = real(transfo.SU(2:end,:)); if output starting OVER the surface
%   SU = real(transfo.SU(2:end,:)); if output starting AT the surface


GS = PIVRes.GS;
% SU=(length(PIVRes.zPIV)-(3087-SU+1)/GS) +1; %length(PIVRes.zPIV) == size(u,1);
SU=(3087-SU+1)/GS; %length(PIVRes.zPIV) == size(u,1);
% (3087-SU+1) is right side up
% (3087-SU+1)/GS is right side up in PIV Water
% length(PIVRes.zPIV)-(3087-SU+1)/GS) +1 is up side down in PIV Water
PIVSurf = PIVRes.PIVW1_Surface ; 

%D = [];
uInt = NaN(size(SU));

for col = 1:size(SU,2)
    
    ggu = u(:,col);
    fvp = find(~isnan(ggu), 1, 'first'); % First Velocity Position: Returns
    % the first element that is not NaN, which is the surface.
    
    if isempty(fvp)
        continue;
    end
    
    ggu = ggu(fvp:end); % The velocity from the surface to the end.
    
    zi = fvp - PIVSurf(col) : size(u,1) - PIVSurf(col); % Initial Data Sites: 
    % The location of velocity measurments in the z-direction. 
    hh = SU(:,col);
    zt = hh(hh < size(u,1) + 50) - PIVSurf(col); % Target data sites. Notice
    % that the SU was NEVER flipped. SU is calculated using the FFT of the flipped surface; We added 50 to be sure that we capture
    % everything.
        
    %D = [D (fvp - PFSurf(col))];% Checks whether the location of first PIV
    % measurement is between 0 and 1. We can also check this 
    %D = [D zt(1)]; That's the average target height measurement (should average to 0.5)
    
    % We only extrapolate on the bottom (near the surface) , not on the top
    % (away from surface).
    ggu_Int1 = interp1(zi,ggu,zt, Type, 'extrap'); % when it is specified by 'extrap', interp1 performs extrapolation for out of range values.
    ggu_Int = ggu_Int1;
    ggu_Int(find(zt>max(zi), 1 ):end) = NaN; %NaN for any element of zt that is outside the interval
    % spanned by zi
    uInt(1:length(ggu_Int),col) = ggu_Int;
    
end

 
%%% Mask for 


% if min(D) < 0
%     display 'There are some measurements below the surface!!!';
% end

end
