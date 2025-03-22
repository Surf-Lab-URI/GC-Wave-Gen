function [IM] = CorrectPIVSURFLensDistortion(img)
% This function corrects the lens distortion of the LFV camera by using the
% results of the PIVSurf lens distortion toolbox. 

%%% MATLAB ***less accurate in this case***
% fc = [6664.04724796681 6657.04923701764];
% cc = [2008.11238871688 1600.21758360679];
% alpha_c = [ 0.00000 ] ;
% kc = [ -0.0779025138414783  0.0781228723983685  0.004457779745842  -0.001965164721823  0.200165359762379 ] ;
% MeanReprojectionError = 0.353213749829286;

% Caltech
% Focal Length:          fc = [ 6703.69524   6689.91145 ] +/- [ 5.64394   5.54363 ]
% Principal point:       cc = [ 2007.10055   1629.99570 ] +/- [ 3.08304   3.73487 ]
% Skew:             alpha_c = [ 0.00000 ] +/- [ 0.00000  ]   => angle of pixel axes = 90.00000 +/- 0.00000 degrees
% Distortion:            kc = [ -0.09152   0.15992   0.00337   -0.00187  0.00000 ] +/- [ 0.00093   0.00546   0.00012   0.00010  0.00000 ]
% Pixel error:          err = [ 0.31945   0.27466 ]

fc = [ 6703.69524   6689.91145 ];
cc = [ 2007.10055   1629.99570 ];
alpha_c = [ 0.00000 ];
kc = [ -0.09152   0.15992   0.00337   -0.00187  0.00000 ];
err = [ 0.31945   0.27466 ];

IM1 = img;

dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
fc_new = dist_amount * fc;
KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1);0 fc_new(2) cc(2) ; 0 0 1];

IM = rect(IM1,eye(3),fc,cc,kc,alpha_c,KK_new);

end



function [Irec] = rect(I,R,f,c,k,alpha,KK_new);


if nargin < 5,
   k = [0;0;0;0;0];
   if nargin < 4,
      c = [0;0];
      if nargin < 3,
         f = [1;1];
         if nargin < 2,
            R = eye(3);
            if nargin < 1,
               error('ERROR: Need an image to rectify');
            end;
         end;
      end;
   end;
end;


if nargin < 7,
   if nargin < 6,
		KK_new = [f(1) 0 c(1);0 f(2) c(2);0 0 1];
   else
   	KK_new = alpha; % the 6th argument is actually KK_new   
   end;
   alpha = 0;
end;



% Note: R is the motion of the points in space
% So: X2 = R*X where X: coord in the old reference frame, X2: coord in the new ref frame.


if ~exist('KK_new'),
   KK_new = [f(1) alpha*f(1) c(1);0 f(2) c(2);0 0 1];
end;


[nr,nc] = size(I);

%Irec = 255*ones(nr,nc);
Irec = NaN(nr,nc);


[mx,my] = meshgrid(1:nc, 1:nr);
px = reshape(mx',nc*nr,1);
py = reshape(my',nc*nr,1);

rays = inv(KK_new)*[(px - 1)';(py - 1)';ones(1,length(px))];


% Rotation: (or affine transformation):

rays2 = R'*rays;

x = [rays2(1,:)./rays2(3,:);rays2(2,:)./rays2(3,:)];


% Add distortion:
xd = apply_distortion(x,k);


% Reconvert in pixels:

px2 = f(1)*(xd(1,:)+alpha*xd(2,:))+c(1);
py2 = f(2)*xd(2,:)+c(2);


% Interpolate between the closest pixels:

px_0 = floor(px2);


py_0 = floor(py2);
py_1 = py_0 + 1;


good_points = find((px_0 >= 0) & (px_0 <= (nc-2)) & (py_0 >= 0) & (py_0 <= (nr-2)));

px2 = px2(good_points);
py2 = py2(good_points);
px_0 = px_0(good_points);
py_0 = py_0(good_points);

alpha_x = px2 - px_0;
alpha_y = py2 - py_0;

a1 = (1 - alpha_y).*(1 - alpha_x);
a2 = (1 - alpha_y).*alpha_x;
a3 = alpha_y .* (1 - alpha_x);
a4 = alpha_y .* alpha_x;

ind_lu = px_0 * nr + py_0 + 1;
ind_ru = (px_0 + 1) * nr + py_0 + 1;
ind_ld = px_0 * nr + (py_0 + 1) + 1;
ind_rd = (px_0 + 1) * nr + (py_0 + 1) + 1;

ind_new = (px(good_points)-1)*nr + py(good_points);



Irec(ind_new) = a1 .* I(ind_lu) + a2 .* I(ind_ru) + a3 .* I(ind_ld) + a4 .* I(ind_rd);



return;


% Convert in indices:

fact = 3;

[XX,YY]= meshgrid(1:nc,1:nr);
[XXi,YYi]= meshgrid(1:1/fact:nc,1:1/fact:nr);

%tic;
Iinterp = interp2(XX,YY,I,XXi,YYi); 
%toc

[nri,nci] = size(Iinterp);


ind_col = round(fact*(f(1)*xd(1,:)+c(1)))+1;
ind_row = round(fact*(f(2)*xd(2,:)+c(2)))+1;

good_points = find((ind_col >=1)&(ind_col<=nci)&(ind_row >=1)& (ind_row <=nri));
end

function [xd,dxddk] = apply_distortion(x,k)


% Complete the distortion vector if you are using the simple distortion model:
length_k = length(k);
if length_k <5 ,
    k = [k ; zeros(5-length_k,1)];
end;


[m,n] = size(x);

% Add distortion:

r2 = x(1,:).^2 + x(2,:).^2;

r4 = r2.^2;

r6 = r2.^3;


% Radial distortion:

cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;

if nargout > 1,
	dcdistdk = [ r2' r4' zeros(n,2) r6'];
end;


xd1 = x .* (ones(2,1)*cdist);

coeff = (reshape([cdist;cdist],2*n,1)*ones(1,3));

if nargout > 1,
	dxd1dk = zeros(2*n,5);
	dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
	dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;
end;


% tangential distortion:

a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta_x = [k(3)*a1 + k(4)*a2 ;
   k(3) * a3 + k(4)*a1];

aa = (2*k(3)*x(2,:)+6*k(4)*x(1,:))'*ones(1,3);
bb = (2*k(3)*x(1,:)+2*k(4)*x(2,:))'*ones(1,3);
cc = (6*k(3)*x(2,:)+2*k(4)*x(1,:))'*ones(1,3);

if nargout > 1,
	ddelta_xdk = zeros(2*n,5);
	ddelta_xdk(1:2:end,3) = a1';
	ddelta_xdk(1:2:end,4) = a2';
	ddelta_xdk(2:2:end,3) = a3';
	ddelta_xdk(2:2:end,4) = a1';
end;

xd = xd1 + delta_x;

if nargout > 1,
	dxddk = dxd1dk + ddelta_xdk ;
    if length_k < 5,
        dxddk = dxddk(:,1:length_k);
    end;
end;


return;

% Test of the Jacobians:

n = 10;

lk = 1;

x = 10*randn(2,n);
k = 0.5*randn(lk,1);

[xd,dxddk] = apply_distortion(x,k);


% Test on k: OK!!

dk = 0.001 * norm(k)*randn(lk,1);
k2 = k + dk;

[x2] = apply_distortion(x,k2);

x_pred = xd + reshape(dxddk * dk,2,n);


norm(x2-xd)/norm(x2 - x_pred)
end

