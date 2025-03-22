function Y = vec_int2(x,Bound,L_poly,L_trans)

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
Length = (max([x1,Bound(1)-L_poly,101]):min([Bound(2)+L_poly,length(x)-100,x2]));
idx = Bound(1):Bound(2);
L2 = Length(~ismember(Length,idx));
SS = x(L2);
fit1 = polyfit(L2(~isnan(SS)),SS(~isnan(SS)),10); %polynomial of 10th order
PolySurf2 = polyval(fit1,Length);
PolySurf3 = x;
PolySurf3(Length) = PolySurf2;
PolySurf3 = PolySurf3 - mean(PolySurf3([Bound(1)-1,Bound(2)+1])) + mean(x([Bound(1)-1,Bound(2)+1]));
a = x(1:Bound(1)-1);
c = x(Bound(2)+1:end);
PS_mod = PolySurf3(Bound(1)-1:Bound(2)+1);
LinIntrp = interp1([Bound(1)-1,Bound(2)+1],[x(Bound(1)-1),x(Bound(2)+1)],Bound(1)-1:Bound(2)+1);
w1L = 1:-1/(L_trans-1):0;
b1 = w1L.*LinIntrp(1:L_trans)+(1-w1L).*PS_mod(1:L_trans);
b2 = (1-w1L).*LinIntrp(end-L_trans+1:end)+w1L.*PS_mod(end-L_trans+1:end);
b= [b1 PS_mod(L_trans+1:end-L_trans) b2];
b([1,end])= [];

Y = [a b c] ;