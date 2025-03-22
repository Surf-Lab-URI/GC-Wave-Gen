% This MATLAB script computes the velocity field by using the two fused PIV
% images. It basically employs the function "Compute Velocities" to perform
% the PIV.


clearvars -except expNameList FI LI; clc;

DataPath = 'D:\';

% expNameList = {'4'};
if ~exist('expNameList', 'var')
    prompt = 'Experiment Name List (e.g. 4,5,6 or 4) = ';
    expNameList = strsplit(input(prompt, 's'), ','); % The number of exper-
    % iments. For a multiple numbers of experiments use comma as delimiter.
end

for expIndex = 1:length(expNameList) % Experiment Loop
    
    expName = expNameList{expIndex};
    
    LoadPathPIV = [DataPath 'Data\Experiment ' expName '\Processed Raw Images\PIV Fused Images\'];
    LoadPathPix = [DataPath 'Data\Experiment ' expName '\Processed Raw Images\PixRes\'];
    
    SavePathRaw = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIV Velocities - Raw\'];
    SavePathPRE = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIVRes\'];
    SavePathRNO = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIV Velocities - Raw No Outliers\'];
    SavePathPPR = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIV Velocities - Post Processed\'];
    
    if ~exist(SavePathRaw, 'dir') || ~exist(SavePathPRE, 'dir') || ... 
       ~exist(SavePathRNO, 'dir') || ~exist(SavePathPPR, 'dir')
        mkdir(SavePathRaw); mkdir(SavePathPRE);
        mkdir(SavePathRNO); mkdir(SavePathPPR);
    end
    
    PixResDir = dir([LoadPathPix '*.mat']);
    
    % FI = 1; LI = length(PIVFusedDir);
    if ~exist('FI', 'var') || ~exist('LI', 'var')
        
        disp(['The length of directory is ' num2str(length(PIVFusedDir)) '.']);
        prompt = 'File Numbers (e.g. all or 10-15 or 5) = ';
        NOF = input(prompt, 's');
        
        if strfind(NOF, 'all')
            FI = 1;
            LI = length(PIVFusedDir);
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
    
    IntrWndw = [128 64 32 16 8]; % Interrogation window sizes
    GrdSpc = [64 32 16 8 4]; % Grid spacing sizes
    
    for index = FI:LI % Main Loop
        
        imNum = PixResDir(index).name(end-7:end-4);
        
        load([LoadPathPIV 'Exp' expName '_PIVInput_' imNum '.mat']);
        load([LoadPathPix 'Exp' expName '_pixRes_' imNum '.mat']);
        
        IM1 = PIV.IM_a;
        IM2 = PIV.IM_b;
        
        CompVel = ComputeVelocities_Quick_NoFilt(IM1, IM2, Mask, IntrWndw, GrdSpc);
        
        
        % Create and Save Raw PIV Output
        PIVRes.xPIV = CompVel.xPIV; % The x coordinates of center of IntrWndws
        PIVRes.zPIV = CompVel.zPIV; % The y coordinates of center of IntrWndws
        PIVRes.PF_Surface = (pixRes.PF_Surface(PIVRes.xPIV) - 1)/CompVel.GS;
        PIVRes.GS = CompVel.GS; % Final grid spacing
        
        FileName = ['Exp' expName '_CompVel_' imNum '.mat'];
        save([SavePathRaw FileName], 'CompVel', 'pixRes', 'PIVRes');
        
        FileName = ['Exp' expName '_PIVRes_' imNum '.mat'];
        save([SavePathPRE FileName], 'PIVRes');
        
        
        % Create and Save Raw PIV Output Without Outliers
        Cartesian = RemoveOutliers_fab(CompVel);
        Cartesian.w = RemoveOutliers(CompVel.delta_z, CompVel.dcor);
        Cartesian.Mask = CompVel.Mask;
        
        FileName = ['Exp' expName '_CartVel_' imNum '.mat'];
        save([SavePathRNO FileName], 'Cartesian', 'pixRes', 'PIVRes');
        
        
        % Create and Save Post Processed (Smoothed) PIV Output
        delx = CompVel.delta_x; % Displacement, or total velocity in the x- 
        % direction. Vectors going downwind (left to right) are positive.
        dely = CompVel.delta_z; % Displacement, or total velocity in the y- 
        % direction. vectors going away from surface (upward) are positive.
        Mask = CompVel.Mask;
        
        delx(CompVel.dcor < 0.4) = NaN;
        dely(CompVel.dcor < 0.4) = NaN;
        
        for i = 1:size(delx,2)
            m = CompVel.Mask(:,i);
            p = find(isnan(m), 1, 'first'); % Finds the first NaN element
            for k = 1:length(m)-p
                
                if (p-k) <= 0
                    delx(p+k-1,i) = delx(p+k-2,i);
                    dely(p+k-1,i) = dely(p+k-2,i); % Fills in with the same
                    % value as above , when wave height is higher than half
                    % the image.
                else
                    delx(p+k-1,i) = -delx(p-k,i);
                    dely(p+k-1,i) = -dely(p-k,i);
                end
                
            end
        end
        
        % Smooth with considering the first 3 and last 2 columns of NaNs.
        delx = smoothn(delx, 0.4, 'robust');
        dely = smoothn(dely, 0.4, 'robust');
        Cartesian.u = delx;
        Cartesian.w = dely;
        Cartesian.Mask = CompVel.Mask;
        
        % Smooth without considering the first 3 and last 2 column of NaNs.
        % delx = smoothn(delx(:,3:end-2), 0.4, 'robust');
        % dely = smoothn(dely(:,3:end-2), 0.4, 'robust');
        % Cartesian.u = [NaN(size(Mask,1),2) delx NaN(size(Mask,1),2)];
        % Cartesian.w = [NaN(size(Mask,1),2) dely NaN(size(Mask,1),2)];
        % Cartesian.Mask = CompVel.Mask;
        
        FileName = ['Exp' expName '_CartVel_' imNum '.mat'];
        save([SavePathPPR FileName], 'Cartesian', 'pixRes', 'PIVRes');
        
        disp(['Pair ' imNum ' done.']);
        
    end
    
end
