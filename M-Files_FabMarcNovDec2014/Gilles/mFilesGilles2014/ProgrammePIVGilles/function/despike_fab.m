function [f]=despike_fab(s);

n=0;
ip1=1;
ip2=[1 1];
f=s;
x_i=[1:length(s)];
x=find((~isnan(f)));
f=f(x);
f_i=interp1(x,f,x_i, 'linear', 'extrap');
f=f_i;


while(length(ip2)~=length(ip1))

ip1=ip2;
[f, ip2] = func_despike_phasespace3d(f);

n=n+1;
x_i=[1:length(s)];
x=find((~isnan(f)));
f=f(x);
f_i=interp1(x,f,x_i, 'linear', 'extrap');
f=f_i;

end
end

