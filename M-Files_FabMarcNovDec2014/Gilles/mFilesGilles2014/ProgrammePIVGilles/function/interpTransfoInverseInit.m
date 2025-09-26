function imgTransfoInverse = interpTransfoInverseInit(surf,pivRes,ImageEntree,SU)
%% Operation inverse de interpTransfo
sPSa = surf.z_s(end:-1:1).*pivRes;   
imPSa = surf.img;
sPS = size(imPSa,1)*pivRes - sPSa;
uui = nan(size(ImageEntree));
SU=SU-mean(SU(1,:))+mean(sPS);
for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
    xxt=[2048-1:-1:0];
    uui(:,col) = interp1(xii,ImageEntree(1:length(xii),col),xxt,'spline',NaN); 
end
imgTransfoInverse=uui;
