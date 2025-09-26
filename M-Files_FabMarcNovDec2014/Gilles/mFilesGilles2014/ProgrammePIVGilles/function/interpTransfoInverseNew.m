function imgTransfo = interpTransfoInverseNew(ImageEntree,SU)

imgTransfo = nan(size(ImageEntree));
for col = 1:size(ImageEntree,2)
    xii = SU(:,col);
    xxt=[size(ImageEntree,1):-1:1];
    imgTransfo(:,col) = interp1(xxt,ImageEntree(1:length(xii),col),xii,'spline',NaN);
end
