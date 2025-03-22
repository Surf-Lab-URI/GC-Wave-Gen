function [row, col, coefs, const, Coef] = coef_RWallV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

x_max = x(end);
tol = 10^(-14);

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
z_num = length(z);

im1 = mod(i-1 - 1,x_num)+1;
ip1 = mod(i-1 + 1,x_num)+1;
if x(im1) < x(i)
    xIm1 = x(im1);
else
    xIm1 = x(im1) - x_max;
end
if x(ip1) > x(i)
    xIp1 = x(ip1);
else
    xIp1 = x(ip1) + x_max;
end

%quadratic interpolation of surface
b = [Surf(im1); Surf(i); Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);


%Compute PhiX
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = Roots <= xIp1 & Roots >= x(i);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root S not found at point (%d, %d)', i, k)
else
    PhiX = (Roots(1) - x(i))/dx;
end

%Find intersection between surface and z-gridlines
xS = x(i)+PhiX*dx;
zS = z(k);

% Find slope at intersection point
slopeS =  (2*coef(1)*(xS)+coef(2));
    
%compute angle of normal relative to x-axis
gamma_s = angleT(1, -1/slopeS);

%find intersection of normal and x-grid line at given point (point A) -----
mN = - 1/slopeS;
bN = zS - mN*xS;

xA = x(i);
zA = mN*xA+bN;
    
%set as default interpolation point (xM, xM)
xOpt = xA;
zOpt = zA;

%Define x relevant for determining distance
w = 1:z_num-k-1;
c1 = z(k + w - 1) <= zOpt & z(k + w) > zOpt;
ind1 = find(c1,1);
clear c1 w

%Define indices of additional point used points used in interpolation
if isempty(ind1) && zOpt > z(end-1)
    ind1 = z_num - (k-1) - 1;
end
if isempty(ind1)
    fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%Cordinates of points used in linear interpolation
% top
iIndT = i;
kIndT = k + ind1;

% bottom
iIndB = i;
kIndB = k + ind1 - 1;

% look for potential diagonal intersection (point C)-----------------------
% number of grid points to look through
searchNum = max(-(floor(-dx/dz)), 1);
if ~periodic
    searchNum = min(searchNum, i);
end

% point type of possible grid points
points = pType(ip1,k:k+searchNum);
ind2 = find(points(2:max(2,length(points)-1))~=-1,1);
if isempty(ind2)
    ind2 = 0;
end
%compare two possible approximations to see which one is better
%if diagonal better, update xOpt to diagonal

if ind2 <= ind1 && ind2 ~= 0
    xT = xIp1;
    zT = z(k + ind2);
    mD = (z(k) - zT)/(x(i) - xT);
    bD = z(k) - mD*x(i);
    
    %find intersection between diagonal and normal
    xCzC = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xC = xCzC(1);
    zC = xCzC(2);
    xOpt = xC;
    zOpt = zC;
    
    %Coordinates of points used in linear interpolation
    % bottom
    iIndB = i;
    kIndB = k;
    % top
    iIndT = ip1;
    kIndT = k + ind2;
end

%Now, construct coefficients
Coef = zeros(x_num, z_num);
Coef(i,k) = -2/dz^2 - 2/(PhiX*dx^2);
Coef(i, k + 1) = 1/dz^2;
Coef(i, k - 1) = 1/dz^2;
Coef(im1 , k) = 2/((PhiX + 1)*dx^2);
Coef(ip1 , k) = 2/(PhiX*(PhiX + 1)*dx^2);

dSOpt = distP([xS, zS], [xOpt, zOpt]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dTOpt = distP([xOpt, zOpt], [x(iIndT), z(kIndT)]);
dBOpt = distP([xOpt, zOpt], [x(iIndB), z(kIndB)]);

%Add interpolated values to appropriate coefficient terms
const = Coef(ip1, k)*dPdnS * dSOpt;
Coef(iIndT, kIndT) = Coef(iIndT, kIndT) + dBOpt/(dBOpt + dTOpt) * Coef(ip1, k);
Coef(iIndB, kIndB) = Coef(iIndB, kIndB) + dTOpt/(dBOpt + dTOpt) * Coef(ip1, k);


%Remove right value
Coef(ip1, k) = 0;

%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

idx = sub2ind(size(Coef),row,col);
coefs = Coef(idx);
