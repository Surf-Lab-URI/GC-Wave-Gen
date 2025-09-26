function imgTransfo = interpTransfo(surf1,pivRes,ImageEntree,SU)
%% Interpole ImageEntree de facon a ce que la surface definie par surf1 devienne la premiere ligne imgTransfo. Une attenuation en prfondeur est obtenue grace a SU
sPSa = surf1.z_s(end:-1:1).*pivRes;
imPSa = surf1.img;
sPS = size(imPSa,1)*pivRes - sPSa;
uui = nan(size(ImageEntree));
SU=SU-mean(SU(1,:))+mean(sPS);
for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
    xxt=[2048-1:-1:0];
    uui(:,col) = interp1(xxt,ImageEntree(1:length(xxt),col),xii,'spline',NaN); 
end
imgTransfo=uui;
