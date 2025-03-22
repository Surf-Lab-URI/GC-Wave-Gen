function [uInt] = inverseTransformVelField_decay_WATER(u, PIVRes, SU)
% This function is vertically transforming zeta to the z, e.g. the velocity
% u(x, zeta) will be transformed to the velocity u(x, z). We need to notice
% that the transformation matrix should be considered  from the second row,
% since the first line is zeta = 0, i.e. the surface exactly. Moreover, the
% u does not need to be flipped upside down.
%
%   u = uInt; (surface at the top when using imagesc)
%   SU = real(transfo.SU(2:end,:)); 


GS = PIVRes.GS;
SU=(length(PIVRes.zPIV)-(3091-SU+1)/GS) +1; %length(PIVRes.zPIV) == size(u,1);
% (3091-SU+1) is right side up
% (3091-SU+1)/GS is right side up in PIV 
% length(PIVRes.zPIV)-(3091-SU+1)/GS) +1 is up side down in PIV
% zPIV = PIVRes.zPIV/GS;
zPIV = length(PIVRes.zPIV)-PIVRes.zPIV/GS +1;
uInt = NaN(size(SU));

for col = 1:size(u,2)
    
    ggu = u(:,col);
    ggu = ggu(~isnan(ggu));
    
    if isempty(ggu)
        continue;
    end
    
    zi = SU(:, col);
    zi = zi(~isnan(u(:, col))); % Initial data sites, which is the location
    % of velocity measurements in the zeta direction.
    
    Diff_Vec = zPIV - PIVRes.PF_Surface(col);
    zt = zPIV(Diff_Vec < 0); % Target data sites.
    
    ggu_Int1 = interp1(zi, ggu, zt, 'spline', NaN);
    ggu_Int2 = interp1(zi, ggu, zt, 'spline', 'extrap');
    ggu_Int2(find(~isnan(ggu_Int1),1, 'last') + 1:end) = NaN; % ggu_Int1(find(~isnan(ggu_Int2),1, 'last') + 1:end);
    ggu_Int2(1:find(~isnan(ggu_Int1), 1, 'first') - 1) = NaN; % ggu_Int1(1:find(~isnan(ggu_Int2), 1, 'first') - 1);
    
    uInt(zt, col) = ggu_Int2;
    
    uTmp = uInt(:, col);
    if ~isnan(uTmp(ceil(PIVRes.PF_Surface(col)) + 1))
        ggu_IntS = interp1(zi, ggu, zt, 'spline', 'extrap');
        ggu_Int2(1:find(~isnan(ggu_Int1), 1, 'first') - 1) = ggu_IntS(1:find(~isnan(ggu_Int1), 1, 'first') - 1);
        % Using spline interpolation for the first row when it is NaN. This
        % is better than the linear one for the first row near the surface.
        uInt(zt, col) = ggu_Int2;      
    end
    if ~isnan(uTmp(size(uTmp,1) - 3))
        ggu_IntL = interp1(zi, ggu, zt, 'linear', 'extrap');
        ggu_Int2(find(~isnan(ggu_Int1),1, 'last') + 1:end) = ggu_IntL(find(~isnan(ggu_Int1),1, 'last') + 1:end);
        % Using linear interpolation  for the last row when it is NaN. This
        % is better than the spline one for the first row near the surface.
        uInt(zt, col) = ggu_Int2;
    end
    
end

uInt = flipud(uInt);

end
