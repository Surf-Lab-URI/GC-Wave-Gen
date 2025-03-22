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
b = [Surf(im1); Surf(i); Surf(ip1)];

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
        ind2 = find(points(2:max(2,length(points)-1))~=-1,1);
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
        [xC, zC] = linsolve([-mD, 1,-mN, 1], [bD bN]);
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