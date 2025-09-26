function Orbitals_interp = interpOrbitals_lc_lfv2(surf1,pivRes,uu,ww,SU,IntrWndw,GrdSpc)

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


SU=SU-mean(SU(1,:))+mean(sPS);

for col = 1:size(uu,2)
%     keyboard
    xii = SU(:,col);
    xxt=[2048-GS*2:-GS:0];    
    uui(:,col) = interp1(xii,uu(1:length(xii),col),xxt,'spline',NaN);
    wwi(:,col) = interp1(xii,ww(1:length(xii),col),xxt,'spline',NaN);
end
Orbitals_interp.u=uui;
Orbitals_interp.w=wwi;