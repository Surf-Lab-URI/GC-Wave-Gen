function [coefIKp1, coefIK, coefIm1K, const] = coef_BRCorner(x, z, Surf, i, k, P_x, P_z)
dx = x(2)-x(1);
dz = z(2)-z(1);

iOrig = i+1;
x_num = length(x);
im1 = mod(i-1 - 1,x_num)+1;
ip1 = mod(i-1 + 1,x_num)+1;

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1), Surf(i), Surf(ip1)];
i = 1+1;

%quadratic interpolation of surface
b = Surf;
A = [x(i - 1).^2, x(i-1), 1;
x(i).^2, x(i), 1;
x(i + 1).^2, x(i + 1), 1];

%coefficients of surface approximation
coef = linsolve(A,b);

%Compute PhiZ
k = k+1;
PhiZ = -(z(k) - Surf(i))./dz;

%Compute PhiX
Roots = roots(coef(0), coef(1), coef(2) - z(k));

Ind = Roots <= x(i + 1) && Roots >= x(i);
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
    gamma_s = angleT(1, -1/(2*coef(1)*(x(i))+coef(2)));
    gamma_n = angleT(1, -1/(2*coef(1)*(x(i) + deltaX*dx)+coef(2)));
    gamma_r = angleT(1, -1/(2*coef(1)*(x(i) + PhiX*dx)+coef(2)));
    
    %set up matrix as in Noye Eq. 13
    C = [0, dz, 0, 0, dz^2/2;
        -dx, 0, dx^2/2, 0, 0;
        cos(gamma_s), sin(gamma_s), 0, cos(gamma_s)*dz*PhiZ,sin(gamma_s)*dz*PhiZ;
        cos(gamma_n), sin(gamma_n), cos(gamma_n)*deltaX*dx,cos(gamma_n)*deltaZ*dz + sin(gamma_n)*deltaX*dx,sin(gamma_n)*deltaZ*dz,    
        cos(gamma_r), sin(gamma_r), cos(gamma_r)*PhiX*dx,sin(gamma_r)*PhiX*dx, 0];
    
end
if cond(C) > 10^7
    fprintf('Poorly conditioned boundary matrix at (%d, %d) \n', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
    CInv = inv(C);
    
    %Compute coefficients of neighboring interior points
    coefIKp1 = CInv(2+1, 0+1) + CInv(4+1, 0+1);
    coefIK = - CInv(2+1, 0+1) -  CInv(2+1, 1+1) - CInv(4+1, 0+1) - CInv(4+1, 1+1);
    coefIm1K = CInv(2+1, 1+1) + CInv(4+1, 1+1);
    
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