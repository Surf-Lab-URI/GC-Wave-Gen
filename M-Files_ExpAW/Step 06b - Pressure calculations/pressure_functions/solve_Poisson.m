function Out = solve_Poisson(f,Surf,x,z,P_x,P_z,curved,periodic,Ptop)
% ---------------------------------------------------------------------- %
% --------------------------- Start function --------------------------- %
% ---------------------------------------------------------------------- %

% Make sure x is in appropriate format for periodic boundary conditions
x = (x - x(1)) + (x(2) - x(1)); % just apply a shift to x(1) value
x_num  = length(x);
z_num = length(z);

% Essential Parameters --------------------------------------------------

% Find step size between grid points
%  (NOTE MUST BE EQUAL IN X AND Z DIRECTIONS!)
h_x = round(mean(diff(x)),14);   % step size
h_z = round(mean(diff(z)),14);

% Using coordinate consistent indexing X inceases by row and Z by column. 
X = repmat(reshape(x,[x_num,1]),[1,z_num]);
Z = repmat(z,[x_num, 1]);

% Remove rows containing only nan from (Surf,X,Z,f,Ptop) <---- Size of these???
nan_rows = find(all(isnan(f),2));
Surf(nan_rows) = [];
x(nan_rows) = [];
f(nan_rows,:) = [];
X(nan_rows,:) = [];
Z(nan_rows,:) = [];
x_num = length(x);

% Round water surface, convert height to index
surf_index = zeros(x_num,1);
for i = 1:x_num
    surf_index(i) = sum(z<Surf(i));
end
surf_min_index = min(surf_index);

% Cut out columns below surface containing only nan
z = z(surf_min_index+1:end);
f = f(:,surf_min_index+1:end);
X = X(:,surf_min_index+1:end);
Z = Z(:,surf_min_index+1:end);
% Shift so min of surf_index is 0 to match new z above
surf_index = surf_index - surf_min_index;
z_num = length(z);

% Check that f is defined everywhere along the water surface, and if not adjust by taking surface to be at index above
for i = 1:x_num
    if isnan(f(i,surf_index(i)+1))
        surf_index(i) = surf_index(i)+1;
        fprintf('Warning Missing f value above surface at (%d, %d) \n',i, surf_index(i) + surf_min_index);
    end
end

% Classify points ------------------------------------------------------
%
% Constuct 2-Dim'l array pType classifying each point according to its type.
% In order to do this, I add one row of 'ghost points' along the boundary
% For this reason, I briefly expand the system to have dimensions
% (x_num + 2, z_num + 2). The ghost points are eliminated at the end of this 
% section
% 
% See classification scheme in document 'Point_Classification.pdf'
% 
% Standard Points
%     -1 ... exterior point
%     0  ... interior point
%     1  ... left wall
%     2  ... ceiling
%     3  ... right wall
%     4  ... bottom wall
%     5  ... top-left corner
%     6  ... top-right corner
%     7  ... bottom right corner
%     8  ... bottom left corner
% 
% Curved Boundary Points
%     21 ... bottom wall
%     22 ... right wall
%     23 ... left wall
%     24 ... bottom right corner 
%     25 ... bottom left corner
%     26 ... bottom right corner of domain 
%     27 ... bottom left corner of domain
% 

% Briefly expand system to include ghost points
x_num2 = x_num + 2;
z_num2 = z_num + 2;

% Set interior point as default
pType = zeros(x_num2,z_num2);
       
% -----------------------------------
% If boundary condtions non-periodic
% -----------------------------------
if ~periodic

    % Set exterior points to -1
    pType(1, :) = -1;
    pType(x_num2, :) = -1;
    for i = 1:x_num
        pType(i+1, 1:surf_index(i)+1) = -1;
    end

        %Classify boundary points along the bottom surface
        surf_adjust = zeros(x_num,1);

        check = nanmax(abs(diff(surf_index))) + 1;
        % number of points above boundary to include in check. Lower numbers
        %make runtime faster at risk of missing boundary points for steep waves
        for i = 2:x_num+1 
            for k = surf_index(i-1)+2:surf_index(i-1) + check + 1
                if pType(i,k) ~= -1
                    n_bds = sum([pType(i,k-1)==-1, pType(i-1,k)== -1,pType(i,k+1)==-1, pType(i+1,k)== -1]);
                    % Classifying by number of exterior neighbors
                    if n_bds == 0
                        %interior points
                        continue
                    elseif n_bds == 1
                        % Wall Boundary Points
                        if pType(i - 1, k) == -1
                            pType(i, k) = 21;
                        elseif pType(i, k + 1) == -1
                            pType(i,k) = 2;
                        elseif pType(i + 1,k) == -1
                            pType(i, k) = 23;
                        elseif pType(i, k - 1) == -1
                            pType(i,k) = 24;
                        end
                        elseif n_bds == 2
                            % Sub-classify order 2 boundary points
                            if pType(i - 1, k) == -1 && pType(i, k + 1) == -1
                                pType(i,k) = 5;
                            elseif pType(i, k + 1) == -1 && pType(i + 1, k) == -1
                                pType(i, k) = 6;
                            elseif pType(i + 1, k) == -1 && pType(i , k - 1) == -1
                                pType(i,k) = 27;
                            elseif pType(i, k - 1) == -1 && pType(i - 1, k) == -1
                                pType(i, k) = 28;
                            else
                                fprintf("Error Unknown boundary type of order two encountered at (%s,%s)",num2str(i - 1), num2str(k - 1))
                                error('Operation aborted') % sys.exit(1) TO BE ADJUSTED
                            end
                            elseif n_bds == 3
                                pType(i,k) = -1; % set point as exterior point
                                % move surface index up by one
                                surf_adjust(i-1) = surf_adjust(i - 1) + 1;
                                %                        print("Warning Boundary point (%s, %s)"  % (i - 1, k - 1) +\
                                %                              " with 3 exterior neighbors removed")
                            else
                                print("You missed a case")
                                error('Operation aborted') % sys.exit(1) TO BE ADJUSTED
                    end
                end
            end
        end

    % Chop pType to eliminate ghost points
    pType = pType(2:end-1, 2:end-1);

    %Adjust water surface index to new level after points removed
    surf_index = surf_index + surf_adjust;


    %left boundary of domain
    pType(1, surf_index(1)+1:z_num) = 1;

    %right boundary
    pType(x_num, surf_index(x_num)+1:z_num) = 3;

    %top boundary
    pType(1:x_num, z_num) = 2;

    %corners
    pType(1, z_num) = 5;
    pType(x_num, z_num) = 6;
    pType(x_num, surf_index(x_num)+1) = 7;
    pType(1, surf_index(1)+1) = 8;
end


%%% ------------------------------
% If boundary conditions periodic
%%% ------------------------------    

if periodic
    
    % Set exterior points to -1
    pType(1, 1:(surf_index(x_num) + 1)) = -1;
    pType(x_num2, 1:(surf_index(1) + 1)) = -1;
    for i = 1:x_num
        pType(i + 1, 1:(surf_index(i) + 1)) = -1;
    end
    %Classify bottom boundary points
    surf_adjust = zeros(x_num,1);   
    
    check = nanmax(abs(diff(surf_index)));
    check = max(check, abs(surf_index(1) - surf_index(end)));
    check = check + 1; % put 10 if using many nan's in the quality check
    
    % number of points above boundary to include in check. Lower numbers
    %make runtime faster at risk of missing boundary points for steep waves    
    
    for i = 2:x_num+1 
            for k = surf_index(i-1)+2:surf_index(i-1) + check + 1
                if pType(i,k) ~= -1
                    n_bds = sum([pType(i,k-1)==-1, pType(i-1,k)== -1,pType(i,k+1)==-1, pType(i+1,k)== -1]);
                    % Classifying by number of exterior neighbors
                    if n_bds == 0
                        %interior points
                        continue
                    elseif n_bds == 1
                        % Wall Boundary Points
                        if pType(i - 1, k) == -1
                            pType(i, k) = 21;
                        elseif pType(i, k + 1) == -1
                            pType(i,k) = 2;
                        elseif pType(i + 1,k) == -1
                            pType(i, k) = 23;
                        elseif pType(i, k - 1) == -1
                            pType(i,k) = 24;
                        end
                        elseif n_bds == 2
                            % Sub-classify order 2 boundary points
                            if pType(i - 1, k) == -1 && pType(i, k + 1) == -1
                                pType(i,k) = 5;
                            elseif pType(i, k + 1) == -1 && pType(i + 1, k) == -1
                                pType(i, k) = 6;
                            elseif pType(i + 1, k) == -1 && pType(i , k - 1) == -1
                                pType(i,k) = 27;
                            elseif pType(i, k - 1) == -1 && pType(i - 1, k) == -1
                                pType(i, k) = 28;
                            else
                                fprintf("Error Unknown boundary type of order two encountered at (%s,%s)",num2str(i - 1), num2str(k - 1))
                                error('Operation aborted') % sys.exit(1) TO BE ADJUSTED
                            end
                            elseif n_bds == 3
                                pType(i,k) = -1; % set point as exterior point
                                % move surface index up by one
                                surf_adjust(i-1) = surf_adjust(i - 1) + 1;
                                %                        print("Warning Boundary point (%s, %s)"  % (i - 1, k - 1) +\
                                %                              " with 3 exterior neighbors removed")
                            else
                                print("You missed a case")
                                error('Operation aborted') % sys.exit(1) TO BE ADJUSTED
                    end
                end
            end
    end
    % Chop pType to eliminate ghost points               
    pType = pType(2:end-1, 2:end-1);
    
    %Adjust water surface index to new level after points removed
    surf_index = surf_index + surf_adjust;
    
    pType(:, z_num) = 2;
end

%Adjust to curved or non-curved pType
if ~curved
    pType(pType >= 0) = mod(pType(pType >= 0),10);
end

% Check for any interior data points with no forcing value

for i = 1:x_num
    for k = 1:z_num
        if pType(i,k) >= 0 && isnan(f(i,k))
            fprintf('Error forcing matrix incomplete. Missing value at address (%s, %s) \n', num2str(i),num2str(k));
            error('Operation aborted')
        end
    end
end

%  % Define Boundary Conditions --------------------------------------------------
% '''
% Boundary condition are defined by 3 matrices ('Dirichlet', 'P_x', 'P_z')
% of size (x_num, z_num).Each element [i,k] of the matrices represents one 
% point on the grid. The Cartesian coordinates of this point are given by 
% X[i,k], Z[i,k].
% 
% Matrices only need to be filled at grid points of the corresponding boundary
% type. For example, 'Dirichlet' only needs to contain values for every boundary
% point with dirichlet conitions. Superflous values will be ignored.
% 
% Neumann boundary conditions are specified by the matrices 'P_x' and 'P_z'. 
% These represent pressure gradients dP/dx and dP/dz, respectively. They must be
% given at all Neumann boundary points where the surface normal has an x-component
% (for P_x) and a z-component (for P_z). On coner points, both must be filled.
% 
% Examples
%    * Dirichlet = 0*X    %Sets P = 0 at all dirichlet boundary points
% 
% 
%    * P_x = 0*X
%      P_z = 0*Z          %Sets dP/dn 0 = at all Neumann boundary points
%   
%     
%    * for i in range(1,x_num)
%         P_z[i,surf_index[i]] += -AIR_DENSITY*Eta_tt[i]
%                        %Adds vertical wave movement to bottom bdry condition
%  
%     
%    * P_x = -AIR_DENSITY*(np.multiply(U,U_x) +np.multiply(W,U_z))
%      P_z = -AIR_DENSITY*(np.multiply(U,W_x) +np.multiply(W,W_z))  
%                        %Creates bdry condition based enitrely off advection
% 
% '''
% % Neumann-Conditions on side boundaries <-------------- On bottom boundary???
PX = 0*X;
PZ = 0*X;
for i = 1:x_num
    PX(i, surf_index(i)+1) = P_x(x(i));
    PZ(i, surf_index(i)+1) = P_z(x(i));
end

% Dirichlet condition on top boundary
Dirichlet = 0*X;
Dirichlet(:, end) = Ptop;
%     
% % Set up system -------------------------------------------------------
%           
% '''
% Using hypthothetical lexographic enumeration of grid points as follows
%     
%     +z
%     |  0   1   2   3   4
%     |  5   6   7   8   9
%     |  10  11  12  13  14
%     |  15  16  17  18  19
%     |  20  21  22  23  24
%     |______________________+x
%     
% Grid point outside the region are skipped and indexed -99.
% 
% Finite difference equations depend on point type as well as type of boundary
% condition. The coefficient matrix is thus set up by a series of if statemetns.
% For each point type, can switch between Dirichlet and Neumann conditions by 
% uncommenting the appropriate code. In order to change boundary condition on
% a subset of points not already identified by a unique pType, one must first 
% create new pType label. 
% 
% '''
% Set up matrix of indices
Ind = -99*ones(x_num,z_num);
ind = 0;
for k = z_num:-1:1
    for i = 1:x_num
        if pType(i,k) ~= -1
            Ind(i,k) = ind;
            ind = ind+1;
        end
    end
end
dh2_x = 1/(h_x^2);            
dh2_z = 1/(h_z^2);
     
% Set up coefficient matrix
row = []; col = []; coef = []; b = [];
for k = z_num:-1:1
        for i = 1:x_num
            
            im1 = mod((i - 2),x_num)+1;
            ip1 = mod(i,x_num)+1;
            
            %interior points
            if pType(i,k) == 0
                row = [row;repmat(Ind(i,k),5,1)];
                col = [col;[Ind(i,k),Ind(ip1,k),Ind(i,k-1),Ind(im1,k),Ind(i,k+1)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z), dh2_x, dh2_z, dh2_x, dh2_z]'];
                b = [b;-f(i,k)];
                
            %exterior points
            elseif pType(i,k) == -1
                continue
            
            %bdry points on left side
            elseif pType(i,k) == 1
                
                %Neumann Condition
                row = [row;repmat(Ind(i,k),4,1)];
                col = [col;[Ind(i,k),Ind(ip1,k),Ind(i,k-1), Ind(i,k+1)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z), 2*dh2_x, dh2_z, dh2_z]'];
                b = [b;-f(i,k) + 2*PX(i, k)/h_x];
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];
%                 '''
            
            % bdry points on top
            elseif pType(i,k) == 2
%                 '''
%                 %Neumann Condition
%                 row = [row;repmat(Ind(i,k),4,1)];
%                 col = [col;[Ind(i,k),Ind(i,k-1),Ind(ip1,k), Ind(im1,k)]']; 
%                 coef = [coef;[-2*(dh2_x + dh2_z),2*dh2_z, dh2_x, dh2_x]'];
%                 b = [b;-f(i,k) -2*PZ(i,k)/h_z];
%                 
%                 '''
                % For Dirichlet Condition
                row = [row;Ind(i,k)];
                col = [col;Ind(i,k)]; 
                coef = [coef;1];
                b = [b;Dirichlet(i,k)];
                
                
                
            % bdry points on right side
            elseif pType(i,k) == 3
                
                %For Neumann Condition
                row = [row;repmat(Ind(i,k),4,1)];
                col = [col;[Ind(i,k),Ind(im1,k),Ind(i,k-1), Ind(i,k+1)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z),2*dh2_x, dh2_z, dh2_z]'];
                b = [b;-f(i,k) + 2*PX(i, k)/h_x];
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];
%                 '''
                
            %bdry points on bottom
            elseif pType(i,k) == 4
                
                %Neumann Condition
                row = [row;repmat(Ind(i,k),4,1)];
                col = [col;[Ind(i,k),Ind(i,k+1),Ind(im1,k), Ind(ip1,k)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z), 2*dh2_z, dh2_x, dh2_x]'];
                b = [b;-f(i,k) + 2*PZ(i, k)/h_z];
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];
%                 '''
            
            % bottom-right corners
            elseif pType(i,k) == 7
                %Neumann Condition
                row = [row;repmat(Ind(i,k),3,1)];
                col = [col;[Ind(i,k),Ind(im1,k),Ind(i,k+1)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z), 2*dh2_x, 2*dh2_z]'];
                b = [b;-f(i,k) + 2*PZ(i, k)/h_z - 2*PX(i, k)/h_x];
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];
%                 '''
                
            % bottom-left corners
            elseif pType(i,k) == 8
                %Neumann Condition
                row = [row;repmat(Ind(i,k),3,1)];
                col = [col;[Ind(i,k),Ind(ip1,k),Ind(i,k+1)]']; 
                coef = [coef;[-2*(dh2_x + dh2_z),2*dh2_x, 2*dh2_z]'];
                b = [b;-f(i,k) + 2*PX(i, k)/h_x + 2*PZ(i, k)/h_z];
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];
%                 '''
            
            % top-left corners
            elseif pType(i,k) == 5
                
                %Dirichlet Condition
                row = [row;Ind(i,k)];
                col = [col;Ind(i,k)]; 
                coef = [coef;1];
                b = [b;Dirichlet(i,k)];
                
%                 '''           
%                 %For Neumann Condition
%                 row = [row;repmat(Ind(i,k),3,1)];
%                 col = [col;[Ind(i,k),Ind(ip1,k),Ind(i,k-1)]']; 
%                 coef = [coef;[-2*(dh2_x + dh2_z), 2*dh2_x, 2*dh2_z]'];
%                 b = [b;-f(i,k) + 2*PX(i,k)/h_x - 2*PZ(i,k)/h_z];
%                 '''
            
            % top-right corners
            elseif pType(i,k) == 6
%                 '''
%                 %Neumann Condition
%                 row = [row;repmat(Ind(i,k),3,1)];
%                 col = [col;[Ind(i,k),Ind(im1,k),Ind(i,k-1)]']; 
%                 coef = [coef;[-2*(dh2_x + dh2_z), 2*dh2_x, 2*dh2_z]'];
%                 b = [b;-f(i,k) - 2*P_x(x(i))/h_x - 2*P_z(x(i))/h_z];
%                 b.append(-f(i,k) - 2*P_x(x(i))/h_x - 2*P_z(x(i))/h_z)
%                 
%                 '''
                %For Dirichlet Condition
                row = [row;Ind(i,k)];
                col = [col;Ind(i,k)]; 
                coef = [coef;1];
                b = [b;Dirichlet(i,k)];                
            
            % left-wall curved
            elseif pType(i,k) == 21
%                 
                try
                    error('1')
                    %Neumann Condition Method Noye
                    [coefIp1K, coefIKp1, coefIKm1, coefIK, const] = ...
                        coef_LWallV2_Noye(x, z, Surf, i, k, P_x, P_z);
                    
                    row = [row;repmat(Ind(i,k),4,1)];
                    col = [col;[Ind(i,k), Ind(ip1,k), Ind(i,k+1), Ind(i,k-1)]'];
                    coef = [coef;[coefIK, coefIp1K, coefIKp1, coefIKm1]'];
                    b = [b;(-f(i,k) - const)];
                    
                catch
                    %Neumann Condition Method Morton
                    [rowApp, colApp, coefs, const, Coef] =...
                        coef_LWallV2_Morton(x,z,Surf,i,k,P_x,P_z,pType,periodic);
                    
                    row = [row;repmat(Ind(i,k),length(rowApp),1)];
                    idx = sub2ind(size(Ind),rowApp, colApp);
                    col = [col;Ind(idx)];
                    coef = [coef;coefs];
                    b = [b;(-f(i,k) - const)];
                end
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];  
% %                 '''
%             
            % right-wall curved
            elseif pType(i,k) == 23
                try
                    error('1')
                    %Neumann Condition Method Noye
                    [coefIm1K, coefIKp1, coefIKm1, coefIK, const] =...
                        coef_RWallV2_Noye(x, z, Surf, i, k, P_x, P_z);
                    
                    row = [row;repmat(Ind(i,k),4,1)];
                    col = [col;[Ind(i,k), Ind(im1,k), Ind(i,k+1), Ind(i,k-1)]'];
                    coef = [coef;[coefIK, coefIm1K, coefIKp1, coefIKm1]'];
                    b = [b;(-f(i,k) - const)];
                catch
                    %Neumann Condition Method Morton
                    [rowApp, colApp, coefs, const, Coef] =...
                        coef_RWallV2_Morton(x,z,Surf,i,k,P_x,P_z,pType,periodic);
                    
                    row = [row;repmat(Ind(i,k),length(rowApp),1)];
                    idx = sub2ind(size(Ind),rowApp, colApp);
                    col = [col;Ind(idx)];
                    coef = [coef;coefs];
                    b = [b;(-f(i,k) - const)];
                end             
                
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];  
%                 '''
%                 
            %bottom-wall curved
            elseif pType(i,k) == 24
                try 
                    %Neumann Condition Method Noye
                    
                    [coefIKp1, coefIp1K, coefIm1K, coefIK, const] =...
                        coef_BWall_Noye(x, z, Surf, i, k, P_x, P_z);
                    
                    row = [row;repmat(Ind(i,k),4,1)];
                    col = [col;[Ind(i,k), Ind(ip1,k), Ind(im1,k), Ind(i,k+1)]'];
                    coef = [coef;[coefIK, coefIp1K, coefIm1K, coefIKp1]'];
                    b = [b;(-f(i,k) - const)];
                    
                catch
                    %Neumann Condition Method Morton
                    [rowApp, colApp, coefs, const, Coef] =...
                        coef_BWallV2_Morton(x,z,Surf,i,k,P_x,P_z,pType,periodic);
                    
                    row = [row;repmat(Ind(i,k),length(rowApp),1)];
                    idx = sub2ind(size(Ind),rowApp, colApp);
                    col = [col;Ind(idx)];
                    coef = [coef;coefs];
                    b = [b;(-f(i,k) - const)];
                end 
            
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];  
%                 '''

            elseif pType(i,k) == 27
                try 
                    error('1')
                    %Neumann Condition Method Noye
                    [coefIKp1, coefIK, coefIm1K, const] =...
                        coef_BRCorner_Noye(x, z, Surf, i, k, P_x, P_z);
                    
                    row = [row;repmat(Ind(i,k),3,1)];
                    col = [col;[Ind(i,k), Ind(im1,k), Ind(i,k+1)]'];
                    coef = [coef;[coefIK, coefIm1K, coefIKp1]'];
                    b = [b;(-f(i,k) - const)];
                    
                catch
                    %Neumann Condition Method Morton
                    [rowApp, colApp, coefs, const, Coef] =...
                        coef_BRCornerV2_Morton(x,z,Surf,i,k,P_x,P_z,pType,periodic);
                    
                    row = [row;repmat(Ind(i,k),length(rowApp),1)];
                    idx = sub2ind(size(Ind),rowApp, colApp);
                    col = [col;Ind(idx)];
                    coef = [coef;coefs];
                    b = [b;(-f(i,k) - const)];
                end 
             
                
%                 '''
%                 %For Dirichlet Condition
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];  
%                 '''
                
                %bottom left corner
            elseif pType(i,k) == 28
                try
                    %Neumann Condition Method Noye
                    error('1')
                    [coefIKp1, coefIK, coefIp1K, const] =...
                        coef_BLCorner_Noye(x, z, Surf, i, k, P_x, P_z);
                    
                    row = [row;repmat(Ind(i,k),3,1)];
                    col = [col;[Ind(i,k), Ind(ip1,k), Ind(i,k+1)]'];
                    coef = [coef;[coefIK, coefIp1K, coefIKp1]'];
                    b = [b;(-f(i,k) - const)];
                                     
                catch
                    %Neumann Condition Method Morton
                    try
                    [rowApp, colApp, coefs, const, Coef] =...
                        coef_BLCornerV2_Morton(x,z,Surf,i,k,P_x,P_z,pType,periodic);
                    catch
                        keyboard
                    end
                    row = [row;repmat(Ind(i,k),length(rowApp),1)];
                    idx = sub2ind(size(Ind),rowApp, colApp);
                    col = [col;Ind(idx)];
                    coef = [coef;coefs];
                    b = [b;(-f(i,k) - const)];
                end
        

%                 '''
%                 % For Dirichlet Condition
%                 row.append(Ind(i,k))
%                 col.append(Ind(i,k))
%                 coef.append(1.)
%                 b.append(Dirichlet(i,k))
%                 row = [row;Ind(i,k)];
%                 col = [col;Ind(i,k)]; 
%                 coef = [coef;1];
%                 b = [b;Dirichlet(i,k)];  
%                 '''
                
            else
                fprintf('Unknown Point Type Discovered while constructing coefficient matrix at (%d, %d) \n', i,k)
                error('Operation aborted')
            end
        end
end
%Construct sparse matrix and convert to csc format
A = sparse(row+1,col+1,coef); %full(A) to see the full sparse matrix
clear row col coef

% Solve the system -----------------------------------------------------
sol = A\b;

% Rewrite solution in matrix form (from the vector format)
P = nan(x_num,z_num);
for i = 1:x_num
    for k = 1:z_num
        if Ind(i,k) ~= -99
            P(i,k) = sol(Ind(i,k)+1);
        end
    end
end

Out.P = P;
Out.z = z;
Out.x = x;
Out.surf = Surf;
Out.surf_index = surf_index;

% solve_Poisson function is over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%