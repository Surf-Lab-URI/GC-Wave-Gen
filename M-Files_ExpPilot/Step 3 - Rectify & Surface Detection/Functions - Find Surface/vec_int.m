function X = vec_int(x,Bound,L_poly,L_trans)

%-----------------------------------------------------------------------%
%
%%% Input %%%
% x = whole vector
% Bound = External boundaries : [LeftPoint RightPoint];
% L_poly = number of points taken at each side of the boundaries (e.g.,
%          500)
% L_trans = transitional interval (e.g., 50);
%
%%% Output %%%
% X = interpolated vector
%
%-----------------------------------------------------------------------%

X = 1:length(x);
x1 = min(X(~isnan(x)));
x2 = max(X(~isnan(x)));
Length = (max([x1,Bound(1)-L_poly,601]):min([Bound(2)+L_poly,length(x)-600,x2]));
idx = Bound(1):Bound(2);
L2 = Length(~ismember(Length,idx));
SS = x(L2);
fit1 = polyfit(L2(~isnan(SS)),SS(~isnan(SS)),10); %polynomial of 10th order
PolySurf2 = polyval(fit1,Length);
MarkInt2(1) = Bound(1);%-25;
MarkInt2(2) = Bound(2);%+25;
X0 = Bound(1)-(L_poly+1);
PS_mod = PolySurf2([MarkInt2(1):MarkInt2(end)]-X0);
a = x(1:MarkInt2(1)-(L_trans+1));
c = x(MarkInt2(end)+(L_trans+1):end);
IndL = MarkInt2(1)-L_trans:MarkInt2(1);
w1L = (MarkInt2(1)-[IndL])/(MarkInt2(1)-(MarkInt2(1)-L_trans));
w2L = (IndL-(MarkInt2(1)-L_trans))/(MarkInt2(1)-(MarkInt2(1)-L_trans));
bL = w1L.*x(IndL)+w2L.*PolySurf2(IndL-X0);
bL(end)= [];
IndR = MarkInt2(end):MarkInt2(end)+L_trans;
w1R = (IndR-(MarkInt2(end)))/(MarkInt2(end)+L_trans-MarkInt2(end));
w2R = (MarkInt2(end)+L_trans-IndR)/(MarkInt2(end)+L_trans-MarkInt2(end));
bR = w1R.*x(IndR)+w2R.*PolySurf2(IndR-X0);
bR(1) = [];

X = [a bL PS_mod bR c];