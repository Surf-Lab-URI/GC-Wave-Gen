function [coefIKp1, coefIp1K, coefIm1K, coefIK, const] = coef_BWall_Noye(x, z, Surf, i, k, P_x, P_z)

iOrig = i;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i-1 - 1,x_num)+1;
ip1 = mod(i-1 + 1,x_num)+1;

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1); Surf(i); Surf(ip1)];
i = 1+1;
tol = 10^(-14);


%quadratic interpolation of surface
b = Surf;
A = [x(i - 1).^2, x(i-1), 1;
x(i).^2, x(i), 1;
x(i + 1).^2, x(i + 1), 1];

%coefficients of surface approximation
coef = linsolve(A,b);

%Compute PhiZ
PhiZ = -(z(k) - Surf(i))./dz;

%Compute deltaX, deltaZ
deltaX = 0.5*sign(Surf(i+1)-Surf(i-1));
xN = x(i) + deltaX*dx;
deltaZ = -(z(k) - polyval(coef,xN))/dz;

%Compute the angles of each normal
if deltaX < -tol
    gamma_s = angleT(-1, 1/(2*coef(1)*(x(i))+coef(2)));
    gamma_n = angleT(-1, 1/(2*coef(1)*(x(i) + deltaX*dx)+coef(2)));
elseif deltaX > tol
    gamma_s = angleT(1, -1/(2*coef(1)*(x(i))+coef(2)));
    gamma_n = angleT(1, -1/(2*coef(1)*(x(i) + deltaX*dx)+coef(2)));
else
    coefIKp1 = 2/((1 - PhiZ)*dz^2);
    coefIp1K = 1/dx^2;
    coefIm1K = 1/dx^2;
    coefIK = -2/dx^2 - 2/((1 - PhiZ)*dz^2);
    const = -2/((1- PhiZ)*dz) * P_z(x(i));
    fprintf('Warning: Boundary Symmetric at Point (%d, %d)', iOrig, k)
end

% Set up coefficient matrix as in Noye
C = [0, dz, 0, 0, dz^2/2;
    dx, 0, dx^2/2, 0, 0;
    -dx, 0, dx^2/2, 0, 0;
    cos(gamma_s), sin(gamma_s), 0, cos(gamma_s)*dz*PhiZ,sin(gamma_s)*dz*PhiZ;
    cos(gamma_n), sin(gamma_n), cos(gamma_n)*deltaX*dx,cos(gamma_n)*deltaZ*dz + sin(gamma_n)*deltaX*dx,sin(gamma_n)*deltaZ*dz];

if cond(C) > 10^9
    printf('Poorly conditioned boundary matrix at (%d, %d) \n', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

CInv = inv(C);

%Compute coefficinets of each grid point
coefIKp1 = CInv(2+1, 0+1) + CInv(4+1, 0+1);
coefIp1K = CInv(2+1, 1+1) + CInv(4+1, 1+1);
coefIm1K = CInv(2+1, 2+1) + CInv(4+1, 2+1);
coefIK =  - CInv(2+1, 0+1) - CInv(2+1, 1+1) - CInv(2+1, 2+1) -...
CInv(4+1, 0+1) - CInv(4+1, 1+1) - CInv(4+1, 2+1);

xS = x(i);
xN = x(i) + deltaX*dx;
zS = z(k) + PhiZ*dz;
zN = z(k) + deltaZ*dz;

const = (CInv(2+1,3+1) + CInv(4+1, 3+1)).*...
(cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS)) +...
(CInv(2+1,4+1) + CInv(4+1, 4+1)).*...
(cos(gamma_n)*P_x(xN) + sin(gamma_n)*P_z(xN));