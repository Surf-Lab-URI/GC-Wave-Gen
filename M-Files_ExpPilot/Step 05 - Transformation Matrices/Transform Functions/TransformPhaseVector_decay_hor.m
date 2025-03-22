function [PHInt] = TransformPhaseVector_decay_hor(PH, PIVRes, SU)
% This function is horizontally transforming the phases from x to Xi,  e.g.
% the phase PH(x,Zeta) will be transformmed to the PH(Xi,Zeta). Note that,
%  [a,b]=ismember(PIVRes.xPIV,PixRes.Combo_Surface_X);
%   PH  = unwrap(PixRes.Combo_Surface_eta_smth_phase(b)); or  PH=Cartesian.Phase;
%   SU = imag(transfo.SU(1,:));


GS = PIVRes.GS;
SU = SU/GS;
xPIV = PIVRes.xPIV/GS;

PHInt = NaN(size(PH));
    
    ggu = PH;
    
    xi = SU + xPIV; % Initial data sites
    xt = xPIV; % Target data sites. These are Xi.
        
    % We only extrapolate on the bottom (near the surface) , not on the top 
    % (away from the surface).
    ggu_Int1 = interp1(xi,ggu,xt, 'linear', 'extrap');
    ggu_Int2 = interp1(xi,ggu,xt, 'linear', NaN); % For the linear methods,
    % interp1 return NaN for any element of xt that is outside the interval
    % spanned by xi. For other methods or when it is specified by 'extrap',
    % interp1 performs extrapolation for out of range values.
    
    ggu_Int1(find(~isnan(ggu_Int2), 1, 'last') + 1:end) = NaN;
    
    PHInt = ggu_Int1;
    
PHInt=wrapToPi(PHInt);

end
