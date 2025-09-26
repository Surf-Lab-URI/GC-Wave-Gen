function Orbitals_interp = interpOrbitals_lc(surf,pivRes,u,w,iws)

sPS = surf;
%%
uu = u(2:end,1:end);
ww = w(2:end,1:end);
% s2s = s2(lim:lim+3784);
uui = nan(size(uu));
wwi = nan(size(ww));

% xii = ;
for col = 1:size(uu,2)
%     keyboard
    surfTemp = sPS((col-1)*iws+iws/2);
    xii = surfTemp:-pivRes*iws:0;
    xxt = 0:pivRes*iws:surfTemp;
    uui(1:length(xii),col) = spline(xii,uu(1:length(xii),col),xxt);
    wwi(1:length(xii),col) = spline(xii,ww(1:length(xii),col),xxt);
end
Orbitals_interp.u=uui;
Orbitals_interp.w=wwi;