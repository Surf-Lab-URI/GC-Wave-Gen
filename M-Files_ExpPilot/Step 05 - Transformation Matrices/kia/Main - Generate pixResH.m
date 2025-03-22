% This MATLAB script creates the pixResH which contains the phase of smooth
% PIV Fused in the Xi and Zeta coordinates.


clear; clc;

DataPath = 'D:\';

% expNameList = {'4'};
if ~exist('expNameList', 'var')
    prompt = 'Experiment Name List (e.g. 4,5,6 or 4) = ';
    expNameList = strsplit(input(prompt, 's'), ','); % The number of exper-
    % iments. For a multiple numbers of experiments use comma as delimiter.
end

for expIndex = 1:length(expNameList) % Experiment Loop
    
    expName = expNameList{expIndex};
    
    LoadPathTrsfFF = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoFFFields\'];
    LoadPathPIVRes = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIVRes\'];
    LoadPathPixRes = [DataPath 'Data\Experiment ' expName '\Processed Raw Images\PixRes\'];
    SavePathPixRes = [DataPath 'Data\Experiment ' expName '\Processed Raw Images\PixRes\'];
    
    PixResDir = dir([LoadPathPixRes '*.mat']);
    
    % FI = 1; LI = length(PixResDir);
    if ~exist('FI', 'var') || ~exist('LI', 'var')
        
        disp(['The length of directory is ' num2str(length(PixResDir)) '.']);
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
        
        imNum = PixResDir(index).name(end-7:end-4);
        
        load([LoadPathPixRes 'Exp' expName '_pixRes_' imNum '.mat']);
        load([LoadPathPIVRes 'Exp' expName '_PIVRes_' imNum '.mat']);
        load([LoadPathTrsfFF 'Exp' expName '_transfoFF_' imNum '.mat']);
        
        % Horizontal Transformation
        u  = unwrap(pixRes.PF_Phase_smth(PIVRes.xPIV));
        SU = imag(transfoFF.SU(2:end,:));
        PhaseVec = TransformVelField_decay_hor_phase(u, PIVRes, SU); 
        
        pixRes.PF_Phase_smth_H = wrapToPi(PhaseVec); % This function  wraps
        % angle in radians to [-pi pi] , while the function wrapTo2Pi wraps 
        % angle in radians to [0 2*pi].
        
        FileName = ['Exp' expName '_pixRes_' imNum '.mat'];
        save([SavePathPixRes FileName], 'pixRes', 'Mask');
        disp(['Pair ' imNum ' done.']);
        
    end

end
