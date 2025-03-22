
function a = angleT(x,y)
if x>0
    a = atan(y/x);
else
    a = pi+atan(y/x);
end
end

function X = distP(x1,x2)
X = sqrt(sum((x1 - x2).^2));
end

function V = exact_Morton(x, z)
A = cos(5*(x -0.5)) + cos(5*z);
end

function V = P_x_Morton(x)
V = -5*sin(5*(x -0.5));
end

function V = P_z_Morton(z)
V = -5*np.sin(5*z);
end

function [row, col, coefs, const, Coef] = coef_BWallV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

i = i+1;
x_max = x(end);
tol = 10^(-14);

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
z_num = length(z);

im1 = mod((i -1 - 1),x_num)+1;
ip1 = mod((i -1 + 1),x_num)+1;
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
b = [Surf(im1), Surf(i), Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);


%Compute PhiZ
k = k+1;
PhiZ = (z(k) - Surf(i))/dz;


slopeS =  (2*coef(1)*(x(i))+coef(2));
if slopeS > 0
    %compute angle of normal relative to x-axis
    gamma_s = angleT(1, -1/slopeS);
    
    xS = x(i);
    zS = z(k) - PhiZ*dz;
    %find intersection of normal and x-grid line at given point (point A) -----
    mN = - 1/slopeS;
    bN = zS - mN*xS;
    
    zA = z(k);
    xA = (zA-bN)/mN;
    
    %find relative index corresponding to this point
    
    
    
    %set as default interpolation point (xM, xM)
    xOpt = xA;
    zOpt = zA;
    
    %Define x relevant for determining distance
    xRel = x;
    xRel((i+1):end) = xRel(i+1:end) - x_max;
    c1 = circshift(fliplr(xRel),i-1);
    c2 = circshift(fliplr(xRel),i);
    ind1 = find(c1(1:end-1)<xOpt & c2(1:end-1)>=xOpt,1);
    clear c1 c2
    
    %Define indices of additional point used points used in interpolation
    if isempty(ind1)
        fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
        Inp = input('Abort operation? Y,N [N]');
        if Inp == 'Y' || 'y'
            error('Operation aborted')
        end
    end
    %Coordinates of points used in linear interpolation
    % left
    iIndL = mod((i -1 - (ind1)),x_num)+1;
    kIndL = k;
    % right
    iIndR = mod((i - 1- ind1 + 1),x_num)+1;
    kIndR = k;
    
    % look for potential diagonal intersection (point C)-----------------------
    % number of grid points to look through
    searchNum = max(-(floor(-dz/dx)), 1);
    if ~periodic
        searchNum = min(searchNum, i-1);
    end
    
    % point type of possible grid points
    if k - 1 >= 1
        points = pType(mod(i-1:-1:i-1-searchNum,x_num)+1, k - 1);
        ind2 = find(points(2:length(points)-1)~=-1,1);
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
        [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
        xOpt = xC;
        zOpt = zC;
        
        %Coordinates of points used in linear interpolation
        % right
        iIndR = i;
        kIndR = k;
        % left
        iIndL = i - ind2;
        kIndL = k - 1;
    end
else
    %compute angle of normal relative to x-axis
    gamma_s = angleT(-1, 1/slopeS);
    
    xS = x(i);
    zS = z(k) - PhiZ*dz;
    %find intersection of normal and x-grid line at given point (point A) -----
    mN = - 1/slopeS;
    bN = zS - mN*xS;
    
    zA = z(k);
    xA = (zA-bN)/mN;
    
    
    %set as default interpolation point (xM, xM)
    xOpt = xA;
    zOpt = zA;
    
    %find relative index corresponding to this point
    %Define x relevant for determining distance
    xRel = x;
    xRel(1:i-1) = xRel(1:i-1) + x_max;
    w = 0:length(xRel)-1;
    c1 = xRel(mod((i-1)+w-1,x_num)+1) <= xOpt & xRel(mod((i-1)+w,x_num)+1)>xOpt;
    ind1 = find(c1,1)-1;
    clear c1
    
    %Define indices of additional point used points used in interpolation
    if isempty(ind1)
        fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
        Inp = input('Abort operation? Y,N [N]');
        if Inp == 'Y' || 'y'
            error('Operation aborted')
        end
    end
    %Cordinates of points used in linear interpolation
    % left
    iIndL = mod((i-1 + (ind1) - 1),x_num)+1;
    kIndL = k;
    % right
    iIndR = mod((i-1 + ind1),x_num)+1;
    kIndR = k;
    
    % look for potential diagonal intersection (point C)-----------------------
    % number of grid points to look through
    searchNum = max(-floor(-dz/dx), 1);
    if ~periodic
        searchNum = min(searchNum, (x_num - 1) - i);
    end
    % point type of possible grid points for diagonal base (point B)
    if k - 1 >= 0
        points = pType(mod(i-1:i-1 + searchNum,x_num)+1, k - 1);
        ind2 = find(points(2:length(points)-1)~=-1,1);
    else
        ind2 = 0;
    end
    %compare two possible approximations to see which one is better
    %if diagonal better, update xOpt to diagonal
    
    if ind2 <= ind1 && ind2 ~= 0
        %compute diagonal
        if x(mod(i-1 + ind2,x_num)+1) > x(i)
            xB = x(mod(i-1 + ind2+1,x_num)+1);
        else
            xB = x(mod(i-1 + ind2,x_num)+1) + x_max;
        end
        zB = z(k - 1);
        %slope and intercept of diagonal
        mD = (z(k) - zB)/(x(i) - xB);
        bD = z(k) - mD*x(i);
        
        %find instersection between diagonal and normal
        [xC, zC] = linsolve([-mD, 1,-mN, 1], [bDl bN]);
        xOpt = xC;
        zOpt = zC;
        
        %Cordinates of points used in linear interpolation
        % right
        iIndL = i;
        kIndL = k;
        % left
        iIndR = mod((i-1 + ind2),x_num)+1;
        kIndR = k - 1;
    end
end

%Now, construct coefficients
Coef = zeros(x_num, z_num);
Coef(i,k) = -2/dx^2 - 2/(PhiZ*dz^2);
Coef(im1, k) = 1/dx^2;
Coef(ip1, k) = 1/dx^2;
Coef(i, k + 1) = 2/((PhiZ + 1)*dz^2);
Coef(i, k - 1) = 2/(PhiZ*(PhiZ + 1)*dz^2);



dSOpt = distP([xS, zS], [xOpt, zOpt]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dROpt = distP([xOpt, zOpt], [x(iIndR), z(kIndR)]);
dLOpt = distP([xOpt, zOpt], [x(iIndL), z(kIndL)]);

%Add interpolated values to appropriate coefficient terms
const = Coef(i, k - 1)*dPdnS * dSOpt;
Coef(iIndL, kIndL) = Coef(iIndL, kIndL) + dROpt/(dROpt + dLOpt) * Coef(i, k - 1);
Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dLOpt/(dROpt + dLOpt) * Coef(i, k - 1);


%Remove bottom value
Coef(i, k - 1) = 0;


%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

coefs = Coef(row, col);

end

function [row, col, coefs, const, Coef] = coef_RWallV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

i = i+1;
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
b = [Surf(im1), Surf(i), Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);


%Compute PhiX
k = k+1;
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = Roots <= xIp1 && Roots >= x(i);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root S not found at point (%d, %d)', i, k)
else
    PhiX = (roots(1) - x(i))/dx;
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
w = 2:z_num-k;
c1 = z(k + w - 1) <= zOpt & z(k + w) > zOpt;
ind1 = find(c1,1);
clear c1

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
ind2 = find(points(2:length(points)-1)~=-1,1);
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
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
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

coefs = Coef(row, col);

end   

function [row, col, coefs, const, Coef] = coef_LWallV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

i = i+1;
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

%quadratic interpolation of surface
b = [Surf(im1), Surf(i), Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);

%Compute PhiX
k = k+1;
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = roots <= x(i) && roots >= xIm1;
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root S not found at point (%d, %d)', i, k)
else
    PhiX = (x(i)-roots(1))/dx;
end

%Find intersection between surface and z-gridlines
xS = x(i)-PhiX*dx;
zS = z(k);

% Find slope at intersection point
slopeS =  (2*coef(1)*(xS)+coef(2));
    
%compute angle of normal relative to x-axis
gamma_s = angleT(-1, 1/slopeS);

%find intersection of normal and x-grid line at given point (point A) -----
mN = - 1/slopeS;
bN = zS - mN*xS;

xA = x(i);
zA = mN*xA+bN;
    
%set as default interpolation point (xM, xM)
xOpt = xA;
zOpt = zA;

%Define x relevant for determining distance
c1 = z(k + w - 1) <= zOpt & z(k + w) > zOpt;
ind1 = find(c1,1);
clear c1

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

% point type of possible grid points
points = pType(im1,k:k+searchNum);
ind2 = find(points(2:length(points)-1)~=-1,1);
if isempty(ind2)
    ind2 = 0;
end
%compare two possible approximations to see which one is better
%if diagonal better, update xOpt to diagonal

if ind2 <= ind1 && ind2 ~= 0
    xT = xIm1;
    zT = z(k + ind2);
    mD = (z(k) - zT)/(x(i) - xT);
    bD = z(k) - mD*x(i);
    
    %find intersection between diagonal and normal
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xOpt = xC;
    zOpt = zC;
    
    %Coordinates of points used in linear interpolation
    % bottom
    iIndB = i;
    kIndB = k;
    % top
    iIndT = im1;
    kIndT = k + ind2;
end

%Now, construct coefficients
Coef = zeros(x_num, z_num);
Coef(i,k) = -2/dz^2 - 2/(PhiX*dx^2);
Coef(i, k + 1) = 1/dz^2;
Coef(i, k - 1) = 1/dz^2;
Coef(ip1 , k) = 2/((PhiX + 1)*dx^2);
Coef(im1 , k) = 2/(PhiX*(PhiX + 1)*dx^2);

dSOpt = distP([xS, zS], [xOpt, zOpt]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dTOpt = distP([xOpt, zOpt], [x(iIndT), z(kIndT)]);
dBOpt = distP([xOpt, zOpt], [x(iIndB), z(kIndB)]);

%Add interpolated values to appropriate coefficient terms
const = Coef(ip1, k)*dPdnS * dSOpt;
Coef(iIndT, kIndT) = Coef(iIndT, kIndT) + dBOpt/(dBOpt + dTOpt) * Coef(im1, k);
Coef(iIndB, kIndB) = Coef(iIndB, kIndB) + dTOpt/(dBOpt + dTOpt) * Coef(im1, k);


%Remove right value
Coef(im1, k) = 0;

%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

coefs = Coef(row, col);

end   

function [row, col, coefs, const, Coef] = coef_BLCornerV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

i = i+1;
x_max = x(end);
tol = 10^(-14);

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
z_num = length(z);

im1 = mod((i -1 - 1),x_num)+1;
ip1 = mod((i -1 + 1),x_num)+1;

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
b = [Surf(im1), Surf(i), Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);

%Compute PhiZ
k = k+1;
PhiZ = (z(k) - Surf(i))/dz;

xS = x(i);
zS = z(k) - PhiZ*dz;
slopeS =  (2*coef(1)*(x(i))+coef(2));
if slopeS > 0
    fprintf('Error: Positive Slope at BL Corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end
%compute angle of normal relative to x-axis
gamma_s = angleT(-1, 1/slopeS);

xS = x(i);
zS = z(k) - PhiZ*dz;

%find intersection of normal and x-grid line at given point (point A) -----
mN = - 1/slopeS;
bN = zS - mN*xS;

zA = z(k);
xA = (zA-bN)/mN;

%set as default interpolation point (xM, xM)
xOptB = xA;
zOptB = zA;

%Define x relevant for determining distance
xRel = x;
xRel(1:i-1) =   xRel(1:i-1) + x_max;
c1 = xRel(mod((i-1)+(1:length(xRel)-1)-1,x_num)+1) <= xOptB & xRel(mod((i-1)+(1:length(xRel)-1),x_num)+1)>zOptB;
ind1 = find(c1,1)-1;
clear c1

%Define indices of additional point used points used in interpolation
if isempty(ind1)
    fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%Coordinates of points used in linear interpolation
% left
iIndL = mod((i-1 + (ind1) - 1),x_num)+1;
kIndL = k;
% right
iIndR = mod((i-1 + ind1),x_num)+1;
kIndR = k;

% look for potential diagonal intersection (point C)
% number of grid points to look through
searchNum = max(-(floor(-dz/dx)), (x_num-1) - (i-1));
if ~periodic
    searchNum = min(searchNum, i-1);
end

% point type of possible grid points
if k - 1 >= 1
    points = pType(mod(i-1:i-1+searchNum,x_num)+1, k - 1);
    ind2 = find(points~=-1,1);
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
    if mod(x(i + ind2),x_num)+1 > x(i)
        xB = x(mod(x(i + ind2),x_num)+1);
    else
        xB = x(mod(x(i + ind2),x_num)+1) + x_max;
    end
    zB = z(k - 1);
    mD = (z(k) - zB)/(x(i) - xB);
    bD = z(k) - mD*x(i);
    
    %find intersection between diagonal and normal
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xOpt = xC;
    zOpt = zC;
    
    %Coordinates of points used in linear interpolation
    % right
    iIndL = i;
    kIndL = k;
    % left
    iIndR = mod(x(i + ind2),x_num)+1;
    kIndR = k - 1;
end

% Compute Necessary Terms for Side Point -------------------

%Compute PhiX
k = k+1;
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = roots <= x(i) && roots >= xIm1;
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root W not found at point (%d, %d)', i, k)
else
    PhiX = (x(i)-roots(1))/dx;
end

%Find intersection between surface and z-gridlines
xW = x(i)-PhiX*dx;
zW = z(k);

% Find slope at intersection point
slopeW =  (2*coef(1)*(xW)+coef(2));
if slopeW>0
    fprintf('Error: Positive slope at BL corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%compute angle of normal relative to x-axis
gamma_w = angleT(-1, 1/slopeS);

%find intersection of normal and nearest interior x-gridline (point D) --------
mN = - 1/slopeW;
bN = zW - mN*xW;

xD = x(i);
zD = mN*xD+bN;
    
%set as default interpolation point (xOpt, xOpt)
xOptL = xD;
zOptL = zD;

%Define x relevant for determining distance
c1 = z(k + w - 1) <= zOpt & z(k + w) > zOpt;
ind1 = find(c1,1);
clear c1

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

% point type of possible grid points
points = pType(im1,k:k+searchNum);
ind2 = find(points(2:length(points)-1)~=-1,1);
if isempty(ind2)
    ind2 = 0;
end
%compare two possible approximations to see which one is better
%if diagonal better, update xOpt to diagonal

if ind2 <= ind1 && ind2 ~= 0
    xT = xIm1;
    zT = z(k + ind2);
    mD = (z(k) - zT)/(x(i) - xT);
    bD = z(k) - mD*x(i);
    
    %find intersection between diagonal and normal
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xOptO = xC;
    zOptL = zC;
    
    %Coordinates of points used in linear interpolation
    % bottom
    iIndB = i;
    kIndB = k;
    % top
    iIndT = im1;
    kIndT = k + ind2;
end

%Now, construct coefficients
Coef = zeros(x_num, z_num);
Coef(i,k) = -2/(PhiZ*dz^2) - 2/(PhiX*dx^2);
Coef(i, k + 1) = 2/((PhiZ + 1)*dz^2);
Coef(i, k - 1) = 2/(PhiZ*(PhiZ + 1)*dz^2);
Coef(ip1, k) = 2/((PhiX + 1)*dx^2);
Coef(im1, k) = 2/(PhiX*(PhiX + 1)*dx^2);


% Add interpolated values to appropriate coefficient terms for bottom points S
dSOptB = distP([xS, zS], [xOptB, zOptB]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dLOptB = distP([xOptB, zOptB], [x(iIndL), z(kIndL)]);
dROptB = distP([xOptB, zOptB], [x(iIndR), z(kIndR)]);

const = Coef(i, k - 1)*dPdnS * dSOptB;
Coef(iIndL, kIndL) = Coef(iIndL, kIndL) + dROptB/(dROptB + dLOptB) * Coef(i, k - 1);
Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dLOptB/(dROptB + dLOptB) * Coef(i, k - 1);

% Add interpolated values to appropriate coefficient terms for side points W
dWOptL = distP([xW, zW], [xOptL, zOptL]);
dPdnW = cos(gamma_w)*P_x(xW) + sin(gamma_w)*P_z(xW);
dTOptL = distP([xOptL, zOptL], [x(iIndT), z(kIndT)]);
dBOptL = distP([xOptL, zOptL], [x(iIndB), z(kIndB)]);

const = const+Coef(im1, k)*dPdnW * dWOptL;
Coef(iIndT, kIndT) = Coef(iIndT, kIndT) + dBOptL/(dBOptL + dTOptL) * Coef(im1, k);
Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dTOptL/(dBOptL + dTOptL) * Coef(im1, k);

%Remove bottom value
Coef(im1, k) = 0;
Coef(i, k - 1) = 0;


%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

coefs = Coef(row, col);
end

function [row, col, coefs, const, Coef] = coef_BRCornerV2_Morton( x, z, Surf, i, k, P_x, P_z, pType, periodic)

i = i+1;
x_max = x(end);
tol = 10^(-14);

dx = x(2)-x(1);
dz = z(2)-z(1);

x_num = length(x);
z_num = length(z);

im1 = mod((i -1 - 1),x_num)+1;
ip1 = mod((i -1 + 1),x_num)+1;

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
b = [Surf(im1), Surf(i), Surf(ip1)];

A = [xIm1^2, xIm1, 1;
    x(i)^2, x(i), 1;
    xIp1^2, xIp1, 1];

coef = linsolve(A,b);

% Compute Compute Necessary Terms for Bottom Point -------------------

%Compute PhiZ
k = k+1;
PhiZ = (z(k) - Surf(i))/dz;

xS = x(i);
zS = z(k) - PhiZ*dz;

slopeS =  (2*coef(1)*(x(i))+coef(2));
if slopeS < 0
    fprintf('Error: Negative Slope at BR Corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
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
ind1 = find(c1(2:end)<xOptB & c2(2:end)>=xOptB,1)-1;
clear c1 c2

%Define indices of additional point used points used in interpolation
if isempty(ind1)
    fprintf('Error at approximating normal derivative at (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
    if Inp == 'Y' || 'y'
        error('Operation aborted')
    end
end

%Coordinates of points used in linear interpolation
% left
iIndL = mod((i-1 + (ind1)),x_num)+1;
kIndL = k;
% right
iIndR = mod((i-1 + ind1 + 1),x_num)+1;
kIndR = k;

% look for potential diagonal intersection (point C)
% number of grid points to look through
searchNum = max(-(floor(-dz/dx)), 1);
if ~periodic
    searchNum = min(searchNum, i-1);
end

% point type of possible grid points
if k - 1 >= 1
    points = pType(mod(i-1:-1:i-1-searchNum-1,x_num)+1, k - 1);
    ind2 = find(points~=-1,1);
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
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xOpt = xC;
    zOpt = zC;
    
    %Coordinates of points used in linear interpolation
    % right
    iIndL = i;
    kIndL = k;
    % left
    iIndR = i-ind2;
    kIndR = k - 1;
end

% Compute Necessary Terms for Side Point -------------------

%Compute PhiX
k = k+1;
Roots = roots([coef(1), coef(2), coef(3) - z(k)]);
Ind = roots <= xIp1  && roots >= x(i);
Roots = Roots(Ind);
if isempty(Roots)
    fprintf('Error: Root W not found at point (%d, %d)', i, k)
else
    PhiX = (roots(1)-x(i))/dx;
end

%Find intersection between surface and z-gridlines
xW = x(i)+PhiX*dx;
zW = z(k);

% Find slope at intersection point
slopeW =  (2*coef(1)*(xW)+coef(2));
if slopeW<0
    fprintf('Error: Negative slope at BR corner at point (%d, %d)', i, k)
    Inp = input('Abort operation? Y,N [N]');
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
c1 = circshift(z,-k);
c2 = circshift(z,-k-1);
ind1 = find(c1(1:end-k-1)<zOptR & c2(1:end-k-1)>=zOptR,1);
clear c1 c2

%Define indices of additional point used points used in interpolation
if isempty(ind1) && zOptR > z(end-1)
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

% point type of possible grid points
points = pType(ip1,k:k+searchNum);
ind2 = find(points(2:length(points)-1)~=-1,1);
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
    [xC, zC] = linsolve([-mD, 1;-mN, 1], [bD; bN]);
    xOptR = xC;
    zOptR = zC;
    
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
Coef(i, k - 1) = 2/(PhiZ*(PhiZ + 1)*dz^2);
Coef(im1, k) = 2/((PhiX + 1)*dx^2);
Coef(ip1, k) = 2/(PhiX*(PhiX + 1)*dx^2);


% Add interpolated values to appropriate coefficient terms for bottom points S
dSOptB = distP([xS, zS], [xOptB, zOptB]);
dPdnS = cos(gamma_s)*P_x(xS) + sin(gamma_s)*P_z(xS);
dLOptB = distP([xOptB, zOptB], [x(iIndL), z(kIndL)]);
dROptB = distP([xOptB, zOptB], [x(iIndR), z(kIndR)]);

const = Coef(i, k - 1)*dPdnS * dSOptB;
Coef(iIndL, kIndL) = Coef(iIndL, kIndL) + dROptB/(dROptB + dLOptB) * Coef(i, k - 1);
Coef(iIndR, kIndR) = Coef(iIndR, kIndR) + dLOptB/(dROptB + dLOptB) * Coef(i, k - 1);

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
Coef(i, k - 1) = 0;


%positions of elements with non-zero coefficients
[row, col] = find(abs(Coef) > tol);

coefs = Coef(row, col);
end
