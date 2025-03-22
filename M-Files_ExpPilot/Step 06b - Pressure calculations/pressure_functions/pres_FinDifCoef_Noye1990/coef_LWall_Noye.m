function [coefIp1K, coefIKp1, coefIKm1, coefIK, const] = coef_LWall(x, z, Surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i-1 - 1,x_num)+1;
ip1 = mod(i-1 + 1,x_num)+1;

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1), Surf(i), Surf(ip1)];
i = 1+1;


%linear interpolation of surface
m = (Surf(i) - Surf(i-1))/(x(i) - x(i-1));
b = Surf(i) - m*x(i);

%Solve for PhiX
k = k+1;
PhiX = (1/m*(z(k) - Surf(i)))/dx;

%Solve for deltaX and deltaY
xN = 1/(m^2 + 1) * (x(i) - b*m + z(k)*m);
zN = m*xN + b;
deltaX = (xN - x(i))/dx;
deltaZ = (zN - z(k))/dz;

%Solve for gamma_s and gamma_n
gamma_s = angleT(-1, 1/m);
gamma_n = angleT(-1, 1/m);

% Set up coefficient matrix as in Noye
C = [dx, 0, dx^2/2, 0, 0;
    0, dz, 0, 0, dz^2/2;
    0, -dz, 0, 0, dz^2/2;
    cos(gamma_s), sin(gamma_s), cos(gamma_s)*dx*PhiX,sin(gamma_s)*dx*PhiX, 0;
    cos(gamma_n), sin(gamma_n), cos(gamma_n)*dx*deltaX, cos(gamma_n)*dz*deltaZ + sin(gamma_n)*dx*deltaX,sin(gamma_n)*deltaZ*dz];

if cond(C) > 100*min(1/dx^2,1/dz^2)
    printf('Poorly conditioned boundary matrix at (%d, %d) \n', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
    CInv = inv(C);
    
    %Compute coefficinets of each grid point
    coefIp1K = CInv(2+1, 0+1) + CInv(4+1, 0+1);
    coefIKp1 = CInv(2+1, 1+1) + CInv(4+1, 1+1);
    coefIKm1 = CInv(2+1, 2+1) + CInv(4+1, 2+1);
    coefIK =  - CInv(2+1, 0+1) - CInv(2+1, 1+1) - CInv(2+1, 2+1) -...
    CInv(4, 0) - CInv(4+1, 1+1) - CInv(4+1, 2+1);
    
    xS = x(i) + PhiX*dx;
    xN = x(i) + deltaX*dx;
    
    const = (CInv(2+1,3+1) + CInv(4+1, 3+1)).*...
    (cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS)) +...
    (CInv(2+1,4+1) + CInv(4+1, 4+1)).*...
    (cos(gamma_n)*P_x(xN) + sin(gamma_n)*P_z(xN));
end