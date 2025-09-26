function imgTransfoInverse = interpTransfoInverse(surf1,pivResMod,ImageEntree,SU,GrdSpc)

%% Operation inverse de interpTransfo
% sPSa = surf1.z_s(end:-1:1).*pivRes;
sPSa = surf1.z_s(end:-1:1)./pivResMod;
imPSa = surf1.img;
sPS = size(imPSa,1) - sPSa;
uui = nan(size(ImageEntree));
% SU=SU-mean(SU(1,:))+mean(sPS);
SU=SU-mean(SU(1,:))+mean(sPS)*pivResMod;
for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
    xxt=[2048-GrdSpc(end):-GrdSpc(end):GrdSpc(end)];
    uui(:,col) = interp1(xii,ImageEntree(1:length(xii),col),xxt,'spline',NaN); 
end
imgTransfoInverse=uui;
