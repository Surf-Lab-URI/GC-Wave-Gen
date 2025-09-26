function imgTransfo = interpTransfoInverse2(ImageEntree,SU)

for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
    xxt=[2048-1:-1:0];
    uui(:,col) = interp1(xxt,ImageEntree(1:length(xxt),col),xii,'spline',NaN); 
end
imgTransfo=uui;
