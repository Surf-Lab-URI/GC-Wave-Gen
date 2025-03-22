function [coefIm1K, coefIKp1, coefIKm1, coefIK, const] = coef_RWallV2_Noye(x, z, Surf, i, k, P_x, P_z)

x_num = length(x);
z_num = length(z);

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

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

coef = solve(A,b);

%Compute PhiX
k = k+1;
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);

Ind = Roots <= x(i+1) && Roots >= x(i);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root S not found at point (%d, %d)' , iOrig, k)
else
    PhiX = (Roots(1) - x(i))/dx;
    xS = x(i) + PhiX*dx;
    zS = z(k);
end

%Compute deltaX, deltaZ
Roots = roots([coef(1), coef(2), coef(3) - (z(k) - 0.5*dz)]);
Ind = Roots <= x(i + 1) && Roots >= x(i);
Roots = Roots(Ind);
if isempty(roots)
    fprintf('Error: Root N not found at point (%d, %d)', iOrig, k)
else
    deltaX = (roots(1) - x(i))/dx;
    deltaZ = -0.5;
    xN = x(i) + dx*deltaX;
    zN = z(k) + dz*deltaZ;
end

gamma_s = angle(1, -1/(2*coef(1)*xS + coef(2)));
gamma_n = angle(1, -1/(2*coef(1)*xN + coef(2)));

% Set up coefficient matrix as in Noye
C = [-dx, 0, dx^2/2, 0, 0;
    0, dz, 0, 0, dz^2/2;
    0, -dz, 0, 0, dz^2/2;
    cos(gamma_s), sin(gamma_s), cos(gamma_s)*dx*PhiX,sin(gamma_s)*dx*PhiX, 0;
    cos(gamma_n), sin(gamma_n), cos(gamma_n)*dx*deltaX,cos(gamma_n)*dz*deltaZ + sin(gamma_n)*dx*deltaX,sin(gamma_n)*deltaZ*dz];

if cond(C) > 10^5
    printf('Poorly conditioned boundary matrix at (%d, %d) \n', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
    CInv = inv(C);
    
    %Compute coefficients of each grid point
    coefIm1K = CInv(2+1, 0+1) + CInv(4+1, 0+1);
    coefIKp1 = CInv(2+1, 1+1) + CInv(4+1, 1+1);
    coefIKm1 = CInv(2+1, 2+1) + CInv(4+1, 2+1);
    coefIK =  - CInv(2+1, 0+1) - CInv(2+1, 1+1) - CInv(2+1, 2+1) -...
    CInv(4+1, 0+1) - CInv(4+1, 1+1) - CInv(4+1, 2+1);
    
    xS = x(i) + PhiX*dx;
    xN = x(i) + deltaX*dx;
    
    const = (CInv(2+1,3+1) + CInv(4+1, 3+1)).*...
    (cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS)) +...
    (CInv(2+1,4+1) + CInv(4+1, 4+1)).*...
    (cos(gamma_n)*P_x(xN) + sin(gamma_n)*P_z(xN));
end