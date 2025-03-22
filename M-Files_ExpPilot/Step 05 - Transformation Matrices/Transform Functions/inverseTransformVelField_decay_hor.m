function uInt = inverseTransformVelField_decay_hor( u, PIVRes, SUH )
% This function is horizontally transforming Xi to the x, e.g. the velocity
% u(Xi,Zeta) will be transformmed to the u(x,Zeta).
%
%   We need to notice that velocity u is already upside down. Moreover, the
%   transformation matrix should be considered from second row since  first
%   row is zeta = 0, i.e. the surface exactly. Finally, we have to consider
%   the horizontally transformed SU which is a function of Xi and Zeta.
%
%   SU = imag(transfo.SU(2:end,:));
%   u  = TransformVelField_decay_hor(uInt, PIVRes, SU);
%   SUH = TransformSUField_decay_hor(imag(transfo.SU(2:end,:)), PIVRes, SU);


GS = PIVRes.GS;
SUH = SUH/GS; % Change the SU from pixel to PIV resolution.
xPIV = PIVRes.xPIV/GS;

uInt = NaN(size(u));

for row = 1:size(SUH,1)
    
    ggu = u(row,:);
    
    xi = xPIV - SUH(row,:);
    xt = xPIV; % Target data sites. These are x's from 1 to 985.
    
    % We only extrapolate on the bottom (near the surface), not on the top,
    % away from surface.    
    ggu_Int1 = interp1(xi,ggu,xt, 'linear', 'extrap');
    ggu_Int2 = interp1(xi,ggu,xt, 'linear', NaN);
    
    ggu_Int1(find(~isnan(ggu_Int2),1,'last') + 1:end) = NaN;
    ggu_Int1(1:find(~isnan(ggu_Int2),1, 'first') - 1) = NaN; 
    
    % We could put the linearly extrapolated values back:
    % ggu_Int1(find(~isnan(ggu_Int2), 1, 'last') + 1:end) = ggu_Int1(find(~isnan(ggu_Int2), 1, 'last') + 1:end);
    % ggu_Int1(1:find(~isnan(ggu_Int2), 1, 'first') - 1) = ggu_Int1(1:find(~isnan(ggu_Int2), 1, 'first') - 1);
    
    % We could compensate a first or last column of NaN:
    % ggu_Int3 = interp1(xi,ggu,xt, 'linear', 'extrap');
    % if sum(isnan(ggu_Int3)) > 0
    %     ggu_Int_Prev = interp1(xi,ggu,xt, 'previous', 'extrap');
    %     ggu_Int_Next = interp1(xi,ggu,xt, 'next', 'extrap');
    %     ggu_Int1(find(~isnan(ggu_Int3), 1, 'last') + 1:end) = ggu_Int_Prev(find(~isnan(ggu_Int3), 1, 'last') + 1:end);
    %     ggu_Int1(1:find(~isnan(ggu_Int3), 1, 'first') - 1) = ggu_Int_Next(1:find(~isnan(ggu_Int3), 1, 'first') - 1);
    % end
    
    uInt(row,:) = ggu_Int1 ;
    
end

end