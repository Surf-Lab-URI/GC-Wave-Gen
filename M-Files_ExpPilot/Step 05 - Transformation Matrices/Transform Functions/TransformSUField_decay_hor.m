function [SUH] = TransformSUField_decay_hor(u, PIVRes, SU)
% This function is horizontally transforming the SU fields from x to the Xi
% i.e. the imaginary or real part of the SU(x,Zeta) will be transformmed to
% the SU(Xi,Zeta). For the inverse horizontal transform function we have to
% use the SU which is in the Xi and Zeta coordinates, not in x and z. Also,
% we need to put back the linear values not the NaN where the interpolation
% does not work, since we are going to use the SUH in the initial data site
% that cannot accept NaNs.
%
%   We have to note that the u here is, in fact, the imaginary  part of the
%   SU(x, Zeta). It is  also necessary that the transformation matrix to be
%   considered from second line, since the first line is zeta = 0, i.e. the
%   surface exactly.
%   
%   u  = imag(transfo.SU(2:end,:));
%   SU = imag(transfo.SU(2:end,:));


GS = PIVRes.GS;
SU = SU/GS; % Change the SU from pixel to PIV resolution.
xPIV = PIVRes.xPIV/GS;

SUH = NaN(size(u));

for row = 1:size(SU,1)
    
    ggu = u(row,:);
    
    xi = SU(row,:) + xPIV; % Initial data sites, which are Xi at constant x
    % of 1 to 985.
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
    
    % ggu_Int1(find(~isnan(ggu_Int2),1,'last') + 1:end) = NaN;    
    % ggu_Int1(1:find(~isnan(ggu_Int2),1, 'first') - 1) = NaN;
    
    ggu_IntL = interp1(xi,ggu,xt, 'linear', 'extrap');
    ggu_Int1(find(~isnan(ggu_Int2), 1, 'last') + 1:end) = ggu_IntL(find(~isnan(ggu_Int2), 1,'last') + 1:end);
    ggu_Int1(1:find(~isnan(ggu_Int2), 1, 'first') - 1) = ggu_IntL(1:find(~isnan(ggu_Int2), 1, 'first') - 1);
    
    SUH(row,:) = ggu_Int1;
    
end

end
