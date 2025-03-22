clear
close all
clc

tic
%-----------INPUT DATA-----------%   % depending 
% Flags
include_viscous_term = true;
curved = true;
periodic = false;
Ptop = 0;

% Define experiment
Exp = input('Number of the experiment (e.g. 4): ','s');
Scene = input('Number of the scene (e.g. 3): ','s');

homePath = 'D:\URI_EXP\';
addpath(genpath([homePath,'M_files\']))
resPath = [homePath,'data\EXPW4\EXPW4_Scene3\RESULTS\'];
cartPath = [resPath,'CALCULATED_FIELDS\Cartesian Fields\Cart Fields\'];
cartDiffPath = [resPath,'CALCULATED_FIELDS\Cartesian Fields\CartDiff Fields\'];
transfoPath = [resPath,'Transformations\'];
LaplPath = [resPath,'CALCULATED_FIELDS\Cartesian Fields\Lapl Fields\'];
savePath = [resPath,'Pressure\'];
PairsCartDiff = dir([cartDiffPath,'*cartDiff*']); % find all the pairs in the selected scene
NN = length(PairsCartDiff);
load([resPath,'ExpW',Exp,'_Scene',Scene,'_Parameters.mat']);

% Define constants and put into 'params'
params.AIR_DENSITY = 1.204;     % air density [kg/m3] at 20°C
params.DVISCOSITY = 1.825e-5;   % dynamic viscosity of air [kg/m*s] at 20°C
params.g = 9.81;                 % gravitational acceleration in [m/s2]
params.WATER_DEPTH = 0.7;       % Mean water depth in [m]
params.WATER_DENSITY = 1000;    % in [kg/m3]
params.SURFACE_TENSION = 0.074; % surface tension in [N/m]
params.TOLERANCE = 10e-14;      % numerical tolerance

%-------------------------------------------------------------------------
% Start main loop through each data file

% Track starting datetime
disp(['Pressure code started on: ',datestr(datetime('now'))]);

% Loop on all the images pairs
In = 1;
NN = 1;
for i = In:NN
    % Load file and variables ----------------------------------
    fname = PairsCartDiff(i).name;
    load([cartDiffPath,fname],'cartDiff')
    load([cartPath,'ExpW',Exp,'_Scene',Scene,'_CartVel_',fname(end-7:end)],'PIVRes','PixRes')
    Surf = PIVRes.PF_Surface*4*CST.DX;
    x = PIVRes.xPIV*CST.DX;
    z = PIVRes.zPIV*CST.DX;
    U_x = flipud(cartDiff.u_x)/CST.DT;
    U_z = flipud(cartDiff.u_z)/CST.DT;
    W_z = flipud(cartDiff.w_z)/CST.DT;
    W_x = flipud(cartDiff.w_x)/CST.DT;
    dx = x(4)-x(3);

    % Compute f from gradients
    rho_a = params.AIR_DENSITY;
    f = rho_a*(-2*U_x.*W_z + 2*U_z.*W_x);
    f = transpose(f);

    % Long surface measurements
    xLong = PixRes.Combo_Surface_X(3:4:end-4)*CST.DX;
    surfLong = (3032-PixRes.Combo_Surface_eta_smth(3:4:end-4))*CST.DX;

    % Define pressure bottom boundary condition ----------------
    % Find surface normal acceleration
    lambda_cut = 0.01;   % cutoff wavelength in [m] (same units as xLong)

    % Use Jeff's version of Christoph's function ---------------
    [utSmth, wtSmth, xU] = surface_accel_Fabio(surfLong,xLong,lambda_cut,params); 
%     [utSmth, wtSmth, xU] = surface_accel(surfLong,xLong,lambda_cut,params); % original Jeff Carpenter function
    
    % Restrict computed accelerations to PIV field of view -----
    % Identify index on surfLong where PIV field begins and ends
    baseU = find(xU>=0, 1);
    endU = baseU + length(x)-1;
    x_accel = xU(baseU:endU);
    utSurf = utSmth(baseU:endU);
    wtSurf = wtSmth(baseU:endU);

    % Check that PIV field is correctly identified
    if min(abs(x - x_accel) - dx > params.TOLERANCE)
        disp('---> Error: Acceleration surface incorrectly identified') 
    else
        disp('Surface acceleration: check')
    end
    
    % Viscous term: -------------------------------------------- 
    %%%%%%%%% Insert Laplacian of u (we can calculate directly here)
    if include_viscous_term 
        % Load viscous term data:
        Lapl = [LaplPath,'ExpW',Exp,'_Scene',Scene,'_Lapl_',fname(end-7:end)];
        if exist(Lapl,'file')
            % If Laplacian already computed, just load it
            load(Lapl);
            u_diff_sum = squeeze(L.u_diff_sum')/CST.DX/CST.DT;
            w_diff_sum = squeeze(L.w_diff_sum')/CST.DX/CST.DT;
        else
            % If Laplacian doesn't exist, compute Laplacian
            load([transfoPath,'ExpW',Exp,'_Scene',Scene,'_transfo_',fname(end-7:end)]);
            dirX = [0;1]; % directional derivative along x
            dirZ = [1;0]; % directional derivative along z
            L.u_xx = csapsDiff1(cartDiff.u_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of u ONLY in x direction
            L.u_zz = csapsDiff1(cartDiff.u_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of u ONLY in z direction
            L.w_xx = csapsDiff1(cartDiff.w_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of w ONLY in x direction
            L.w_zz = csapsDiff1(cartDiff.w_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of w ONLY in z direction
            L.u_zz = -L.u_zz; % change sign because z is flipped
            L.w_zz = -L.w_zz; % change sign because z is flipped
            L.Lapl_u = L.u_xx + L.u_zz;
            L.Lapl_w = L.w_xx + L.w_zz;
            L.Lapl_uInt = TransformVelField_decay_0( flipud(L.Lapl_u.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) );
            L.Lapl_wInt = TransformVelField_decay_0( flipud(L.Lapl_w.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) );
            
            % Use the Laplacian extrapolated exactly on the surface
            L.u_diff_sum = L.Lapl_uInt(1,:);
            L.w_diff_sum = L.Lapl_wInt(1,:);
            L.DX = CST.DX;
            L.DT = CST.DT;
            L.x_pix = transfo.x_pix;
            L.zeta_pix = transfo.zeta_pix;
            clear transfo
            
            if ~exist(LaplPath,'dir')
                mkdir(LaplPath)
            end
            save([LaplPath,'ExpW',Exp,'_Scene',Scene,'_Lapl_',fname(end-7:end)],'L')
            
            u_diff_sum = squeeze(L.u_diff_sum')/CST.DX/CST.DT;
            w_diff_sum = squeeze(L.w_diff_sum')/CST.DX/CST.DT;
        end
        disp('Viscous boundary term included')
    else
        u_diff_sum = 0;
        w_diff_sum = 0;
        disp('Viscous boundary term not included')
    end

    %  Pressure gradient components needed for bottom BC --------
    mu = params.DVISCOSITY;
    rho_a = params.AIR_DENSITY;
    P_x = @(X) interp1(x_accel,-rho_a*utSurf + mu*u_diff_sum,X,'spline','extrap');
    P_z = @(Z) interp1(x_accel,-rho_a*wtSurf + mu*w_diff_sum,Z,'spline','extrap');

    % Create pressure field ------------------------------------
    Out = solve_Poisson(f,Surf,x,z,P_x,P_z,curved,periodic,Ptop);
    % unpack the output:
    p2 = Out.P;       % pressure field [Pa]
    % fields: 'P','z','x','surf','surf_index'
    clear Out

    % Return p to original dimensions
    p = nan(length(x),length(z));
    pz_num = size(p2,2);
    p(:,length(z) - pz_num+1:length(z)) = p2;
    clear p2 pznum

    % Save as mat file -------------------------
    press.p = p;
    press.x = x;
    press.z = z;
    filenum = [];
    fname = [savePath,'ExpW',Exp,'_Scene',Scene,'_pressure_',fname(end-7:end)];
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    save(fname,'press','-mat')
    clear cartDiff PIVRes L PixRes
    
    disp(['file finished: ',num2str(i),' of ',num2str(length(PairsCartDiff))])
end

disp(['Pressure code completed on: ',datestr(datetime('now'))]);