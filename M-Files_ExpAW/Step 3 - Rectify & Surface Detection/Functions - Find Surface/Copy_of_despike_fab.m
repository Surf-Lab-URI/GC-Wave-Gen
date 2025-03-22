function [f]=Copy_of_despike_fab(s,opt)

ii=isnan(s);

f=s;
x_i=[1:length(s)];
x=find((~isnan(f)));
f=f(x);
f_i=interp1(x,f,x_i, 'nearest', 'extrap');
f=f_i;



ip1=1;
ip2=[1 1];
while(length(ip2)~=length(ip1))
ip1=ip2;
[f, ip2] = func_despike_phasespace3d(f,0,2); 
%n=n+1;
x_i=[1:length(f)];
x=find((~isnan(f)));
f=f(x);
f_i=interp1(x,f,x_i, 'nearest', 'extrap');
f=f_i;
end

f=f';

if opt==0
    f(ii)=NaN;
end







function [fo, ip] = func_despike_phasespace3d( fi, i_plot, i_opt );
%======================================================================
%
% Version 1.00
%
% This subroutine excludes spike noise from Acoustic Doppler 
% Velocimetry (ADV) data using phasce-space method by 
% Modified Goring and Nikora (2002) method by Nobuhito Mori (2005).
% 
%======================================================================
%
% Input
%   fi     : input data with dimension (n,1)
%   i_plot : =9 plot results (optional)
%   i_opt : = 0 or not specified  ; return spike noise as NaN
%           = 1            ; remove spike noise and variable becomes shorter than input length
%           = 2            ; interpolate NaN using cubic polynomial
%
% Output
%   fo     : output (filterd) data
%   ip     : excluded array element number in fi
%
% Example: 
%   [fo, ip] = func_despike_phasespace3d( fi, 9 );
%     or
%   [fo, ip] = func_despike_phasespace3d( fi, 9, 2 );
%
%
%======================================================================
% Terms:
%
%       Distributed under the terms of the terms of the BSD License
%
% Copyright:
%
%       Nobuhito Mori
%           Disaster Prevention Research Institue
%           Kyoto University
%           mori@oceanwave.jp
%
%========================================================================
%
% Update:
%       1.11    2009/06/09 Interpolation has been added.
%       1.01    2009/06/09 Minor bug fixed
%       1.00    2005/01/12 Nobuhito Mori
%
%========================================================================

nvar = nargin;
if nvar==1
    i_opt  = 0;
    i_plot = 0;
elseif nvar==2
    i_opt = 0;
end

%
% --- initial setup
%

% number of maximum iternation
n_iter = 20;
n_out  = 999;

n      = size(fi,1);
f_mean = nanmean(fi);
f      = fi - f_mean;
lambda = sqrt(2*log(n));

if nargin==1
  i_plot = 0;
end

%
% --- loop
%

n_loop = 1;

while (n_out~=0) & (n_loop <= n_iter)

%
% --- main
%

% step 0
f = f - nanmean(f);
%nanstd(f)

% step 1: first and second derivatives
%ff=filtfilt(ones(1,5)/5,1,f);
f_t  = gradient(f); %f_t=filtfilt(ones(1,5)/5,1,f_t);
f_tt = gradient(f_t); %f_tt=filtfilt(ones(1,5)/5,1,f_tt);
% step 2: estimate angle between f and f_tt axis
if n_loop==1
  theta = atan2( sum(f.*f_tt), sum(f.^2) );
end

% step 3: checking outlier in the 3D phase space
[xp,yp,zp,ip,coef] = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta);

%
% --- excluding data
%

n_nan_1 = size(find(isnan(f)==1),1);
f(ip)  = NaN;
n_nan_2 = size(find(isnan(f)==1),1);
n_out   = n_nan_2 - n_nan_1;

%
% --- end of loop
%

n_loop = n_loop + 1;

end

%
% --- post process
%

go = f + f_mean;
ip = find(isnan(go));

if n_loop < n_iter
  
else
  disp(...
    ['>> Number of outlier   = ', num2str(size(find(isnan(f)==1),1)), ...
     ' : Number of iteration = ', num2str(n_loop-1), ' !!! exceed maximum value !!!'] ...
  )
end

%
% --- interpolation or shorten NaN data
%

if abs(i_opt) >= 1
	% remove NaN from data
    inan = find(~isnan(go));
    fo = go(inan);
    % interpolate NaN data
    if abs(i_opt) == 2
        x   = find(~isnan(go));
        y   = go(x);
        xi  = 1:max(length(fi));
        fo = interp1(x, y, xi, 'linear')';
    end
else
    % output despiked value as NaN
    fo = go;
end

%
% --- for check and  plot
%

if i_plot == 9 

%theta/pi*180
F    = fi - f_mean;
F_t  = gradient(F); 
F_tt = gradient(F_t);
RF = [ cos(theta) 0  sin(theta); 0 1 0 ; -sin(theta) 0 cos(theta)];
RB = [ cos(theta) 0 -sin(theta); 0 1 0 ;  sin(theta) 0 cos(theta)];

% making ellipsoid data
a = coef(1);
b = coef(2);
c = coef(3);
ne  = 32;
dt  = 2*pi/ne;
dp  = pi/ne;
t   = 0:dt:2*pi;
p   = 0:dp:pi;
n_t = max(size(t));
n_p = max(size(p));

% making ellipsoid
for it = 1:n_t
  for is = 1:n_p
    xe(n_p*(it-1)+is) = a*sin(p(is))*cos(t(it));
    ye(n_p*(it-1)+is) = b*sin(p(is))*sin(t(it));
    ze(n_p*(it-1)+is) = c*cos(p(is));
  end
end
xer = xe*RB(1,1) + ye*RB(1,2) + ze*RB(1,3);
yer = xe*RB(2,1) + ye*RB(2,2) + ze*RB(2,3);
zer = xe*RB(3,1) + ye*RB(3,2) + ze*RB(3,3);

% plot figures
figure(1);clf
plot3(f,f_t,f_tt,'b*','MarkerSize',3)
hold on
  plot3(F(ip),F_t(ip),F_tt(ip),'ro','MarkerFaceColor','r','MarkerSize',5)
  plot3(xer,yer,zer,'k-');
hold off
axis equal
grid on
xlabel('u');
ylabel('\Delta u');
zlabel('\Delta^2 u');


figure(2);clf
plot(fi,'k-');
hold on
plot(ip,fi(ip),'ro');
if i_opt==2
    plot(fo,'r-');
end
hold off
%pause

end

function [xp,yp,zp,ip,coef] = func_excludeoutlier_ellipsoid3d(xi,yi,zi,theta);
%======================================================================
%
% Version 1.01
%
% This program excludes the points outside of ellipsoid in two-
% dimensional domain
%
% Input
%   xi : input x data
%   yi : input y data
%   zi : input z data
%   theta  : angle between xi and zi
%
% Output
%   xp : excluded x data
%   yp : excluded y data
%   zp : excluded y data
%   ip : excluded array element number in xi and yi
%   coef : coefficients for ellipsoid
%
% Example: 
%   [xp,yp,zp,ip,coef] = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta);
%
%
%======================================================================
% Terms:
%
%       Distributed under the terms of the terms of the BSD License
%
% Copyright: 
%
%       Nobuhito Mori, Kyoto University
%
%========================================================================
%
% Update:
%       1.01    2009/06/09 Nobuhito Mori
%       1.00    2005/01/12 Nobuhito Mori
%
%========================================================================

%
% --- initial setup
%

n = max(size(xi));
lambda = sqrt(2*log(n));

xp = [];
yp = [];
zp = [];
ip = [];

%
% --- rotate data
%

%theta = atan2( sum(xi.*zi), sum(xi.^2) );

if theta == 0
  X = xi;
  Y = yi;
  Z = zi;
else
  R = [ cos(theta) 0  sin(theta); 0 1 0 ; -sin(theta) 0 cos(theta)];
  X = xi*R(1,1) + yi*R(1,2) + zi*R(1,3);
  Y = xi*R(2,1) + yi*R(2,2) + zi*R(2,3);
  Z = xi*R(3,1) + yi*R(3,2) + zi*R(3,3);
end

%test
%plot3(xi,yi,zi,'b*')
%hold on
%  plot3(X,Y,Z,'r*')
%hold off
%pause

%
% --- preprocess
%

a = lambda*nanstd(X);
b = lambda*nanstd(Y);
c = lambda*nanstd(Z);

%
% --- main
%

m = 0;
for i=1:n
  x1 = X(i);
  y1 = Y(i);
  z1 = Z(i);
  % point on the ellipsoid
  x2 = a*b*c*x1/sqrt((a*c*y1)^2+b^2*(c^2*x1^2+a^2*z1^2));
  y2 = a*b*c*y1/sqrt((a*c*y1)^2+b^2*(c^2*x1^2+a^2*z1^2));
  zt = c^2* ( 1 - (x2/a)^2 - (y2/b)^2 );
  if z1 < 0
    z2 = -sqrt(zt);
  elseif z1 > 0
    z2 = sqrt(zt);
  else
    z2 = 0;
  end

  % check outlier from ellipsoid
  dis = (x2^2+y2^2+z2^2) - (x1^2+y1^2+z1^2);
  if dis < 0 
    m = m + 1;
    ip(m) = i;
    xp(m) = xi(i);
    yp(m) = yi(i);
    zp(m) = zi(i);
  end
end

coef(1) = a;
coef(2) = b;
coef(3) = c;








