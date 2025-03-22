function [coefIKp1, coefIK, coefIp1K, const] = coef_BLCorner(x, z, surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i-1 - 1,x_num)+1;
ip1 = mod(i-1 + 1,x_num)+1;

x = [x(i) - dx, x(i), x(i) + dx];
surf = [surf(im1), surf(i), surf(ip1)];
i = 1+1;

%quadratic interpolation of surface
b = surf;
A = [x(i - 1).^2, x(i-1), 1;
    x(i).^2, x(i), 1;
    x(i + 1).^2, x(i + 1), 1];

coef = solve(A,b);

%Compute PhiZ
k = k+1;
PhiZ = -(z(k) - surf(i))/dz;

%Compute PhiX
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);

Ind = Roots <= x(i) && Roots >= x(i-1);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root R not found at point (%d, %d)' , iOrig, k)
else
    PhiX = (Roots(1) - x(i))/dx;
    
    deltaX = 0.5*PhiX;
    xN = x(i) + deltaX*dx;
    deltaZ = -(z(k) - polyval(coef,xN))/dz;
    %
    
    %compute the angles of each normal for each point
    gamma_s = angleT(-1, 1/(2*coef(1)*(x(i))+coef(2)));
    gamma_n = angleT(-1, 1/(2*coef(1)*(x(i) + deltaX*dx)+coef(2)));
    gamma_r = angleT(-1, 1/(2*coef(1)*(x(i) + PhiX*dx)+coef(2)));    
    
    C = [0, dz, 0, 0, dz^2/2;
        dx, 0, dx^2/2, 0, 0;
        cos(gamma_s), sin(gamma_s), 0, cos(gamma_s)*dz*PhiZ,sin(gamma_s)*dz*PhiZ;
        cos(gamma_n), sin(gamma_n), cos(gamma_n)*deltaX*dx,cos(gamma_n)*deltaZ*dz + sin(gamma_n)*deltaX*dx,sin(gamma_n)*deltaZ*dz,    
        cos(gamma_r), sin(gamma_r), cos(gamma_r)*PhiX*dx,sin(gamma_r)*PhiX*dx, 0]; 
end
if cond(C) > 10^9
    fprintf('Poorly conditioned boundary matrix at (%d, %d)', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
    CInv = inv(C);
    
    %Compute coefficients of neighboring interior points
    coefIKp1 = CInv(2+1, 0+1) + CInv(4+1, 0+1);
    coefIK = - CInv(2+1, 0+1) -  CInv(2+1, 1+1) - CInv(4+1, 0+1) - CInv(4+1, 1+1);
    coefIp1K = CInv(2+1, 1+1) + CInv(4+1, 1+1);
    
    xS = x(i);
    xN = x(i) + deltaX*dx;
    xR = x(i) + PhiX*dx;
    
    %compute constant
    const = (CInv(2+1, 2+1) + CInv(4+1,2+1)).*...
        (cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS)) + ...
        (CInv(2+1,3+1) + CInv(4+1, 3+1)).*...
        (cos(gamma_n)*P_x(xN) + sin(gamma_n)*P_z(xN)) +...
        (CInv(2+1,4+1) + CInv(4+1, 4+1)).*...
        (cos(gamma_r)*P_x(xR) + sin(gamma_r)*P_z(xR));
end