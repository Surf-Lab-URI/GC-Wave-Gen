function uInt = transformVelField_decay_forFab( u, pivRes, SU )
%
% inspired from E:\MFiles\transformation\transformVelField18.m
%
z_PIV = pivRes.zPIV;
%
uInt = nan(size(u));
% SU  = real(transfo.SU);
%
for col = 1:size(SU,2)
    %     col = 500
    mask=pivRes.mask(:,col);
    fvp = find(~isnan(mask),1,'first'); % first velocity position
    if isempty(fvp)
        continue;
    end
    ggu = u(:,col);
    ggu = ggu(fvp:end);
    zi  = z_PIV(fvp):pivRes.GS:z_PIV(end); % initial data sites
    hh  = SU(:,col);
    zt  = hh(hh < max(z_PIV)+150);  % target sites
    %% only extrapolate on the bottom (near the surface), not on the top (away from surface)
    ggu_int1 = interp1(zi,ggu,zt,'linear', 'extrap');
    ggu_int2 = interp1(zi,ggu,zt,'linear', nan);
    ggu_int = ggu_int1;
    ggu_int(find(~isnan(ggu_int2), 1, 'last')+1:end) = nan;
    %%
    uInt(1:length(ggu_int),col) = ggu_int;
end
%
end

