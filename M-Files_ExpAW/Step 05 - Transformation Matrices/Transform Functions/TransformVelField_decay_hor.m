function [uInt] = TransformVelField_decay_hor(u, PIVRes, SU)
% This function is horizontally transforming x to the Xi, e.g. the velocity
% u(x,Zeta) will be transformmed to the u(Xi,Zeta). It is inspired from the
% following file: MFiles\transformation\transformVelField18.m.
%
%   We need to note that velocity u is the output of the TransformVelField_
%   decay(u,pivRes,SU), and therefore, it is already upside down. Thus, the
%   u that is fed here is already transformed vertically from z to Zeta. It
%   is also required that the transformation matrix considered from  second
%   line, since the first line is zeta = 0, i.e. the surface exactly.
%   But if we want uINT at the surface, we should consider the first row.
% 
%   u  = TransformVelField_decay(u, pivRes, SU)
%   SU = real(transfo.SU(2:end,:)); if output starting OVER the surface
%   SU = real(transfo.SU(2:end,:)); if output starting AT the surface


GS = PIVRes.GS;
SU = SU/GS; % Change the SU from pixel to PIV resolution.
xPIV = PIVRes.xPIV/GS;

uInt = NaN(size(u));

for row = 1:size(u,1)
    
    ggu = u(row,:);
    
    xi = SU(row,:) + xPIV; % Initial data sites which are Xi at constant x.
    xt = xPIV; % Target data sites, which are Xi. Notice that these are xi,
    % we just want them to go from 1 to 985 (PIVres).
        
    % We only extrapolate on the bottom (near the surface) , not on the top 
    % (away from the surface).
    if sum(isnan(ggu)) > 0
        ggu_Int1 = interp1(xi,ggu,xt, 'linear', 'extrap');
    else
        ggu_Int1 = interp1(xi,ggu,xt, 'spline', 'extrap');
    end
    
    ggu_Int2 = interp1(xi,ggu,xt, 'linear', NaN); % For the linear methods,
    % interp1 return NaN for any element of xt that is outside the interval
    % spanned by xi. For other methods or when it is specified by 'extrap',
    % interp1 performs extrapolation for out of range values.
    
    ggu_Int1(find(~isnan(ggu_Int2),1,'last') + 1:end) = NaN;    
    ggu_Int1(1:find(~isnan(ggu_Int2),1, 'first') - 1) = NaN; 
    
    % We could put linearly extrapolated values back:
    % ggu_IntL = interp1(xi,ggu,xt, 'linear', 'extrap');
    % ggu_Int1(find(~isnan(ggu_Int2), 1, 'last') + 1:end) = ggu_IntL(find(~isnan(ggu_Int2), 1,'last') + 1:end);
    % ggu_Int1(1:find(~isnan(ggu_Int2), 1, 'first') - 1) = ggu_IntL(1:find(~isnan(ggu_Int2), 1, 'first') - 1);
    
    uInt(row,:) = ggu_Int1;
    
end

end