function Orbitals_interp = interpOrbitals_lc_lfv(surf1,pivRes,uu,ww,IntrWndw,GrdSpc)

sPSa = surf1.z_s(end:-1:1).*pivRes;
imPSa = surf1.img;
sPS = size(imPSa,1)*pivRes - sPSa;

GS=GrdSpc(end);
IW=IntrWndw(end);

%%
%uu = u(2:end,floor(lim/iws):floor((lim+3784)/iws));
%ww = w(2:end,floor(lim/iws):floor((lim+3784)/iws));

% s2s = s2(lim:lim+3784);
uui = nan(size(uu));
wwi = nan(size(ww));
col_number=IW/2:GS:(length(sPS)-IW/2);

for col = 1:size(uu,2)
%     keyboard
    surfTemp = sPS(col_number(col)); 
    xii = surfTemp:-GS:0;
    xxt = 0:GS:surfTemp;
    
    uui(1:length(xxt),col) = spline(xii,uu(1:length(xii),col),xxt);
    wwi(1:length(xxt),col) = spline(xii,ww(1:length(xii),col),xxt);
end
Orbitals_interp.u=uui;
Orbitals_interp.w=wwi;