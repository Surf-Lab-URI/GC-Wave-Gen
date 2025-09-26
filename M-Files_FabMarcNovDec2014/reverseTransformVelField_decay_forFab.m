function uInt = reverseTransformVelField_decay_forFab( u, pivRes, SU )
%
z_PIV = pivRes.zPIV;
uInt = nan(size(u));
%
for col = 1:size(u,2)
    mask=pivRes.mask(:,col);
    fvp = find(~isnan(mask),1,'first'); % first velocity position
    if isempty(fvp)
        continue;
    end
    zt  = z_PIV(fvp):pivRes.GS:z_PIV(end); % target data sites
    zi1 = SU(:,col); %initial positions, too long for now
    ggu1 = u(:,col);
    zi = zi1(~isnan(ggu1));
    ggu = ggu1(~isnan(ggu1));
        
    if isempty(ggu)
        continue;
    end
    %     if length(ggu)<length(zt)-1
    %         keyboard;
    %     end
    ggu_int1 = interp1(zi,ggu,zt,'linear', 'extrap');
    ggu_int2 = interp1(zi,ggu,zt,'linear', nan);
    ggu_int = ggu_int1;
    ggu_int(find(~isnan(ggu_int2), 1, 'last')+2:end) = nan;
    uInt(length(ggu_int):-1:1,col) = ggu_int ; %(...-1)/4 to switch from pixRes coordinates to pivRes positions
    
end
uInt=flipud(uInt);
end

