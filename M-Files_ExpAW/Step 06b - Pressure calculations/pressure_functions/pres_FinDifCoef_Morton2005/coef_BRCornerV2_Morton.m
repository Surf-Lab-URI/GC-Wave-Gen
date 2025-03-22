function [row, col, coefs, const, Coef] = coef_BRCornerV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

x_max = x(end);
tol = 10^(-14);

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
z_num = length(z);

im1 = mod((i-1 - 1),x_num)+1;
ip1 = mod((i-1 + 1),x_num)+1;

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

%quadratic interpolation of surface -------------------
b = [Surf(im1); Surf(i); Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);

% Compute Compute Necessary Terms for Bottom Point -------------------

%Compute PhiZ
PhiZ = (z(k) - Surf(i))/dz;

xS = x(i);
zS = z(k) - PhiZ*dz;

slopeS =  (2*coef(1)*(x(i))+coef(2));
if slopeS < 0
    fprintf('Error: Negative Slope at BR Corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]','s');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end
%compute angle of normal relative to x-axis
gamma_s = angleT(1, -1/slopeS);

xS = x(i);
zS = z(k) - PhiZ*dz;

%find intersection of normal and x-grid line at given point (point A) -----
mN = - 1/slopeS;  %%%% check if it's correct
bN = zS - mN*xS;

zA = z(k);
xA = (zA-bN)/mN;

%set as default interpolation point (xM, xM)
xOptB = xA;
zOptB = zA;

%Define x relevant for determining distance
xRel = x;
xRel(i+1:end) =  xRel(i+1:end) - x_max;
c1 = circshift(fliplr(xRel),i);
c2 = circshift(fliplr(xRel),i+1);
ind1 = find(c1(2:end)<xOptB & c2(2:end)>=xOptB,1);
clear c1 c2

%Define indices of additional point used points used in interpolation
if isempty(ind1)
    fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]','s');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%Coordinates of points used in linear interpolation
% left
iIndL = mod((i-1 - (ind1)),x_num)+1;
kIndL = k;
% right
iIndR = mod((i-1 - ind1 + 1),x_num)+1;
kIndR = k;

% look for potential diagonal intersection (point C)
% number of grid points to look through
searchNum = max(-(floor(-dz/dx)), 1);
if ~periodic
    searchNum = min(searchNum, i-1);
end

% point type of possible grid points
if k - 1 >= 1
    points = pType(mod(i-1:-1:i-1-searchNum,x_num)+1, k - 1);
    ind2 = find(points(2:max(2,length(points)-1))~=-1,1);
    if isempty(ind2)
        ind2 = 0;
    end
else
    ind2 = 0;
end

%compare two possible approximations to see which one is better
%if diagonal better, update xOpt to diagonal

if ind2 <= ind1 && ind2 ~= 0
    %compute diagonal
    if x(i - ind2) < x(i)
        xB = x(i - ind2);
    else
        xB = x(i - ind2) - x_max;
    end
    zB = z(k - 1);
    mD = (z(k) - zB)/(x(i) - xB);
    bD = z(k) - mD*x(i);
    
    %find intersection between diagonal and normal
    [xCzC,~] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xC = xCzC(1);
    zC = xCzC(2);
    xOptB = xC;
    zOptB = zC;
    
    %Coordinates of points used in linear interpolation
    % right
    iIndR = i;
    kIndR = k;
    % left
    iIndL = i-ind2;
    kIndL = k - 1;
end

% Compute Necessary Terms for Side Point -------------------

%Compute PhiX
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = Roots <= xIp1  & Roots >= x(i);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root W not found at point (%d, %d)', i, k)
else
    PhiX = (Roots(1)-x(i))/dx;
end

%Find intersection between surface and z-gridlines
xW = x(i)+PhiX*dx;
zW = z(k);

% Find slope at intersection point
slopeW =  (2*coef(1)*(xW)+coef(2));
if slopeW<0
    fprintf('Error: Negative slope at BR corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]','s');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%compute angle of normal relative to x-axis
gamma_w = angleT(1, -1/slopeS);

%find intersection of normal and nearest interior x-gridline (point D) --------
mN = - 1/slopeW;
bN = zW - mN*xW;

xD = x(i);
zD = mN*xD+bN;
    
%set as default interpolation point (xOpt, xOpt)
xOptR = xD;
zOptR = zD;

%Define x relevant for determining distance
c1 = circshift(z,-k+1);
c2 = circshift(z,-k);
ind1 = find(c1(1:end-k-1)<=zOptR & c2(1:end-k-1)>zOptR,1);
clear c1 c2

%Define indices of additional point used points used in interpolation
if isempty(ind1) && zOptR > z(end-1)
    ind1 = z_num - (k-1) - 1;
end
if isempty(ind1)
    fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]','s');
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
    [xCzC,~] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xC = xCzC(1);
    zC = xCzC(2);
    xOptR = xC;
    zOptR = zC;
    clear xCzC
    
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
Coef(i,k) = -2/(PhiZ*dz^2) - 2/(PhiX*dx^2);
Coef(i, k + 1) = 2/((PhiZ + 1)*dz^2);
if k == 1
    Coef(i, end) = 2/(PhiZ*(PhiZ + 1)*dz^2);
else
    Coef(i, k - 1) = 2/(PhiZ*(PhiZ + 1)*dz^2);
end
Coef(im1, k) = 2/((PhiX + 1)*dx^2);
Coef(ip1, k) = 2/(PhiX*(PhiX + 1)*dx^2);


% Add interpolated values to appropriate coefficient terms for bottom points S
dSOptB = distP([xS, zS], [xOptB, zOptB]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dLOptB = distP([xOptB, zOptB], [x(iIndL), z(kIndL)]);
dROptB = distP([xOptB, zOptB], [x(iIndR), z(kIndR)]);

if k == 1
    const = Coef(i, end)*dPdnS * dSOptB;
    Coef(iIndL, kIndL) = Coef(iIndL, kIndL) + dROptB/(dROptB + dLOptB) * Coef(i, end);
    Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dLOptB/(dROptB + dLOptB) * Coef(i, end);
else
    const = Coef(i, k - 1)*dPdnS * dSOptB;
    Coef(iIndL, kIndL) = Coef(iIndL, kIndL) + dROptB/(dROptB + dLOptB) * Coef(i, k - 1);
    Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dLOptB/(dROptB + dLOptB) * Coef(i, k - 1);
end

% Add interpolated values to appropriate coefficient terms for side points W
dWOptR = distP([xW, zW], [xOptR, zOptR]);
dPdnW = cos(gamma_w)*P_x(xW) + sin(gamma_w)*P_z(xW);
dTOptR = distP([xOptR, zOptR], [x(iIndT), z(kIndT)]);
dBOptR = distP([xOptR, zOptR], [x(iIndB), z(kIndB)]);

const = const+Coef(ip1, k)*dPdnW * dWOptR;
Coef(iIndT, kIndT) = Coef(iIndT, kIndT) + dBOptR/(dBOptR + dTOptR) * Coef(ip1, k);
Coef(iIndB, kIndB) = Coef(iIndB, kIndB) + dTOptR/(dBOptR + dTOptR) * Coef(ip1, k);

%Remove bottom value
Coef(ip1, k) = 0;
if k == 1
    Coef(i, end) = 0;
else
    Coef(i, k - 1) = 0;
end


%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

idx = sub2ind(size(Coef),row,col);
coefs = Coef(idx);