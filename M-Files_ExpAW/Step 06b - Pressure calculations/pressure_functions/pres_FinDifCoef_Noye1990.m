% Bunch of functions for coefficients calculation based on Noye1990

function a = angleT(x,y)
    if x>0
        a = atan(y/x);
    else
        a = pi+atan(y/x);
    end
end

function [coefIKp1, coefIK, coefIm1K, const] = coef_BRCorner(x, z, Surf, i, k, P_x, P_z)
dx = x(2)-x(1);
dz = z(2)-z(1);

iOrig = i+1;
x_num = length(x);
im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

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
end

function [coefIKp1, coefIK, coefIp1K, const] = coef_BLCorner(x, z, Surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1), Surf(i), Surf(ip1)];
i = 1+1;

%quadratic interpolation of surface
b = Surf;
A = [x(i - 1).^2, x(i-1), 1;
    x(i).^2, x(i), 1;
    x(i + 1).^2, x(i + 1), 1];

coef = solve(A,b);

%Compute PhiZ
k = k+1;
PhiZ = -(z(k) - Surf(i))/dz;

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
end

function [coefIm1K, coefIKp1, coefIKm1, coefIK, const] = coef_RWall(x, z, Surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1), Surf(i), Surf(ip1)];
i = 1+1;


%linear interpolation of surface
m = (Surf(i+1) - Surf(i))/(x(i+1) - x(i));
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
gamma_s = angleT(1, -1/m);
gamma_n = angleT(1, -1/m);

% Set up coefficient matrix as in Noye
C = [-dx, 0, dx^2/2, 0, 0;
    0, dz, 0, 0, dz^2/2;
    0, -dz, 0, 0, dz^2/2;
    cos(gamma_s), sin(gamma_s), cos(gamma_s)*dx*PhiX,sin(gamma_s)*dx*PhiX, 0;
    cos(gamma_n), sin(gamma_n), cos(gamma_n)*dx*deltaX,cos(gamma_n)*dz*deltaZ + sin(gamma_n)*dx*deltaX,sin(gamma_n)*deltaZ*dz];

if cond(C) > 10^10
    printf('Poorly conditioned boundary matrix at (%d, %d) \n', iOrig, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
    CInv = inv(C);
    
    %Compute coefficinets of each grid point
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
end

function [coefIp1K, coefIKp1, coefIKm1, coefIK, const] = coef_LWall(x, z, Surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

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
    CInv(4+1, 0+1) - CInv(4+1, 1+1) - CInv(4+1, 2+1);
    
    xS = x(i) + PhiX*dx;
    xN = x(i) + deltaX*dx;
    
    const = (CInv(2+1,3+1) + CInv(4+1, 3+1)).*...
    (cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS)) +...
    (CInv(2+1,4+1) + CInv(4+1, 4+1)).*...
    (cos(gamma_n)*P_x(xN) + sin(gamma_n)*P_z(xN));
end
end

function [coefIKp1, coefIp1K, coefIm1K, coefIK, const] = coef_BWall(x, z, Surf, i, k, P_x, P_z)

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

x = [x(i) - dx, x(i), x(i) + dx];
Surf = [Surf(im1), Surf(i), Surf(ip1)];
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
k = k+1;
PhiZ = -(z(k) - Surf(i))./dz;

%Compute deltaX, deltaZ
deltaX = 0.5*sign(Surf(x+1)-Surf(x-1));
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
end
end

function [coefIm1K, coefIKp1, coefIKm1, coefIK, const] = coef_RWallV2(x, z, Surf, i, k, P_x, P_z)

x_num = length(x);
z_num = length(z);

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

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
end

function [coefIp1K, coefIKp1, coefIKm1, coefIK, const] = coef_LWallV2(x, z, Surf, i, k, P_x, P_z)

x_num = length(x);
z_num = length(z);

iOrig = i+1;

dx = x(2)-x(1);
dz = z(2)-z(1);

im1 = mod(i - 1,x_num);
ip1 = mod(i + 1,x_num);

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

Ind = Roots <= x(i) && Roots >= x(i-1);
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
Ind = Roots <= x(i) && Roots >= x(i-1);
Roots = Roots(Ind);
if isempty(roots)
    fprintf('Error: Root N not found at point (%d, %d)', iOrig, k)
else
    deltaX = (roots(1) - x(i))/dx;
    deltaZ = -0.5;
    xN = x(i) + dx*deltaX;
    zN = z(k) + dz*deltaZ;
end

gamma_s = angle(-1, 1/(2*coef(1)*xS + coef(2)));
gamma_n = angle(-1, 1/(2*coef(1)*xN + coef(2)));

% Set up coefficient matrix as in Noye
C = [dx, 0, dx^2/2, 0, 0;
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
    coefIp1K = CInv(2+1, 0+1) + CInv(4+1, 0+1);
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
end