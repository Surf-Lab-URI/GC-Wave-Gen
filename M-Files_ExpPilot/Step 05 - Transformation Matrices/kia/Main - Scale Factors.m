
clearvars -except expNameList FI LI; clc;

DataPath = 'D:\';

% expNameList = {'4'};
if ~exist('expNameList', 'var')
    prompt = 'Experiment Name List (e.g. 4,5,6 or 4) = ';
    expNameList = strsplit(input(prompt, 's'), ','); % The number of exper-
    % iments. For a multiple numbers of experiments use comma as delimiter.
end

for i = 1:length(expNameList) % Experiment Loop
    
    expName = expNameList{i};
    
    LoadPathPIVRes = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIVRes\'];
    LoadPathTransfo = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoFields\'];
    LoadPathTransfoH = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoHFields\'];
    SavePathScaleFac = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\Scale Factors\'];
    
    if ~exist(SavePathScaleFac, 'dir')
        mkdir(SavePathSU);
    end
    
    PIVDir = dir([LoadPathPIVRes '*.mat']);
    
     % FI = 1; LI = length(PixResDir);
    if ~exist('FI', 'var') || ~exist('LI', 'var')
        
        disp(['The length of directory is ' num2str(length(PIVDir)) '.']);
        prompt = 'File Numbers (e.g. all or 10-15 or 5) = ';
        NOF = input(prompt, 's');
        
        if strfind(NOF, 'all')
            FI = 1;
            LI = length(PixResDir);
        else
            NOF = strsplit(NOF, '-');            
            FI = str2double(NOF{1});
            try
                LI = str2double(NOF{2});
            catch
                LI = FI;
            end
        end
        
    end
    
    
    for index = FI:LI % Main Loop
        
        PairNumber = PIVDir(index).name(end-7:end-4);   
        
        load([LoadPathPIVRes 'Exp4_PIVRes_' PairNumber '.mat']);
        load([LoadPathTransfo 'Exp4_transfo_' PairNumber '.mat']);
        load([LoadPathTransfoH 'Exp4_transfoH_' PairNumber '.mat']);
        
        SUR = real(transfo.SU(2:end, :)); % Real part of the transformation
        % matrix, which contains z while rows are constant Zeta and columns
        % are constant x. That means z(x,zeta).
        SUI = imag(transfo.SU(2:end,:)); % Imaginary part of transformation
        % matrix, which contains xi(x,zeta) - x.
        
        SUHR = real(transfoH.SU(1:end, :)); % Real  part  of transformation
        % matrix, which contains z while rows are constant Zeta and columns
        % are constant Xi. That means z(xi,zeta)  
        
        % Calculating dZ(Xi,Zeta)/dXi and dZ(Xi,Zeta)/dZeta
        Z = smoothn(SUHR, 0, 'robust'); % z = z(xi,zeta)
        Z(~isnan(SUHR)) = SUHR(~isnan(SUHR));
        
        [dZ_dXi, dZ_dZeta] = csapsDiff(Z, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        dZ_dXi(isnan(SUHR)) = SUHR(isnan(SUHR));
        dZ_dZeta(isnan(SUHR)) = SUHR(isnan(SUHR));
        
        % Calculating dXi(Xi,Zeta)/dx and dXi(Xi,Zeta)/dz 
        X = repmat(PIVRes.xPIV,509,1);
        Xi = SUI + X; % This matrix contains Xi(x,Zeta) , rows are constant
        % Zeta and columns are constant x.
        
        Xi_XZ = inverseTransformVelField_decay(Xi, PIVRes, real(transfo.SU(2:end,:))); % Xi(x,z)
        
        Xi_XZ_Smth = smoothn(Xi_XZ, 0, 'robust'); % Xi = Xi(x,z)
        Xi_XZ_Smth(~isnan(Xi_XZ)) = Xi_XZ(~isnan(Xi_XZ));
        
        [dXi_dx, dXi_dz] = csapsDiff(Xi_XZ_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV); % dXi(x,z)/dx, dXi(x,z)/dz
        
        dXi_dx(isnan(Xi_XZ)) = Xi_XZ(isnan(Xi_XZ));
        dXi_dz(isnan(Xi_XZ)) = Xi_XZ(isnan(Xi_XZ));
        
        % Transform dXi(x,z)/dx, dXi(x,z)/dz to dXi(Xi,Zeta)/dx, dXi(Xi,Zeta)/dz 
        dXi_dX = TransformVelField_decay(dXi_dx, PIVRes, real(transfo.SU(2:end,:)));
        dXi_dX = TransformVelField_decay_hor(dXi_dX, PIVRes, imag(transfo.SU(2:end,:)));
        dXi_dZ = TransformVelField_decay(dXi_dz, PIVRes, real(transfo.SU(2:end,:)));
        dXi_dZ = TransformVelField_decay_hor(dXi_dZ, PIVRes, imag(transfo.SU(2:end,:)));
        
        % Calculating dX(Xi,Zeta)/dXi and dX(Xi,Zeta)/dZeta 
        X_SF = TransformVelField_decay(X, PIVRes, real(transfo.SU(2:end,:))); % X(x,Zeta)
        X_SFH = TransformVelField_decay_hor(X_SF, PIVRes, imag(transfo.SU(2:end,:))); % X(Xi,Zeta)
        
        X_SFH_Smth = smoothn(X_SFH, 0, 'robust');
        X_SFH_Smth(~isnan(X_SFH)) = X_SFH(~isnan(X_SFH));
        
        [dX_dXi, dX_dZeta] = csapsDiff(X_SFH_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        dX_dXi(isnan(X_SFH)) = X_SFH(isnan(X_SFH)); % not as good as 1./dXi_dX
        dX_dZeta(isnan(X_SFH)) = X_SFH(isnan(X_SFH)); % dx/dzeta = - dz/dxi
        
        % To double check, we have to notice that dX_dXi should be equal to
        % the inverse of dXi_dX. 
        
        %% Test 
        % Following calculations is another check for (Xi,Zeta) coordinates
        % Since (Xi,Zeta) are orthogonal then dXi/dXi = 1 and dXi/dZeta = 0
        
        % Xi_SFH = TransformVelField_decay_hor(Xi, PIVRes, imag(transfo.SU(2:end,:)));
        % This matrix contains Xi(Xi,Zeta) while rows are constant Zeta and
        % the columns of matrix are constant Xi.
        
        % Xi_SFH_Smth = smoothn(Xi_SFH, 0, 'robust'); % Xi = Xi(Xi,Zeta)
        % Xi_SFH_Smth(~isnan(Xi_SFH)) = Xi_SFH(~isnan(Xi_SFH));
        
        % [dXi_dXi, dXi_dZeta] = csapsDiff(Xi_SFH_Smth, 0.001, PIVRes.xPIV, PIVRes.zPIV);
        
        % dXi_dXi(isnan(Xi_SFH)) = Xi_SFH(isnan(Xi_SFH));
        % dXi_dZeta(isnan(Xi_SFH)) = Xi_SFH(isnan(Xi_SFH));
        
        %% Scale Factors
        
        dX_dXi = 1./dXi_dX;
        dZ_dXi = -dX_dZeta; % Cauchy-Riemann condition dx/dzeta = - dz/dxi,
        % but dX_dZeta is interpolating better at the edges.
        
        h1 = sqrt( (dX_dXi).^2 + (dZ_dXi).^2 );
        h3 = sqrt( (dX_dZeta).^2 + (dZ_dZeta).^2 );
        
        J = (dX_dXi .* dZ_dZeta) - (dX_dZeta .* dZ_dXi); % Jacobian(Xi,Zeta)
        % Jacobian for (Xi, Zeta) coordinates, in theory, it is defined as,
        % J = (dX/dXi)(dZ/dZeta) - (dX/dZeta)(dz/dXi). The second  term  is 
        % small and can be neglected. Note that if we just use the (X,Zeta)
        % coordinates then the Jacobian is J = dZ/dZeta, or better yet , we
        % can cacluate Jacbian analystically at the same time as transfo.SU
        % (i.e. use ak calculations but with ALL modes), which is,
        % J = 1 - akcos(kx)exp(-kZeta).
        
        %% Saving
        FileName = ['Exp' expName '_ScaleFactors_' PairNumber '.mat'];
        save([SavePathScaleFac FileName], 'h1', 'h3', 'J');
        
        disp(['Pair ' PairNumber ' is done.']);
        
    end
        
end
