function imgTransfo255472 = interpTransfo255472(surf1,ImageEntree,SU,GrdSpc)

sPSa = surf1.z_s(end:-1:1);
imPSa = surf1.img;
sPS=size(imPSa,1) - sPSa;
uui = nan(size(ImageEntree));
SU=SU-mean(SU(1,:))+mean(sPS);

for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
%     xxt=[2048-GrdSpc(end):-GrdSpc(end):GrdSpc(end)];
xxt=[2048-GrdSpc(end)*2:-GrdSpc(end):0];
    uui(:,col) = interp1(xxt,ImageEntree(1:length(xxt),col),xii,'spline',NaN); 
end
imgTransfo255472=uui;
