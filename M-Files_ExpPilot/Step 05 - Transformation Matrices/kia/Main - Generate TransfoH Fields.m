

clearvars -except expNameList FI LI; clc; % warning('off', 'all');

DataPath = 'D:\';

% expNameList = {'4'};
if ~exist('expNameList', 'var')
    prompt = 'Experiment Name List (e.g. 4,5,6 or 4) = ';
    expNameList = strsplit(input(prompt, 's'), ','); % The number of exper-
    % iments. For a multiple numbers of experiments use comma as delimiter.
end

for expIndex = 1:length(expNameList) % Experiment Loop
    
    expName = expNameList{expIndex};
    
    LoadPathPIVRes = [DataPath 'Data\Experiment ' expName '\PIV Velocities\PIVRes\'];
    LoadPathTransfo = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoFields\'];
    SavePathTransfoH = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoHFields\'];
    LoadPathTransfoFF = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoFFFields\'];    
    SavePathTransfoFFH = [DataPath 'Data\Experiment ' expName '\Transformation Matrices\TransfoFFHFields\'];
    
    if ~exist(SavePathTransfoH, 'dir') || ~exist(SavePathTransfoFFH, 'dir')
        mkdir(SavePathTransfoH); mkdir(SavePathTransfoFFH);
    end
    
    TransfoDir = dir([LoadPathTransfo '*.mat']);
    
    % FI = 1; LI = length(TransfoDir);
    if ~exist('FI', 'var') || ~exist('LI', 'var')
        
        disp(['The length of directory is ' num2str(length(TransfoDir)) '.']);
        prompt = 'File Numbers (e.g. all or 10-15 or 5) = ';
        NOF = input(prompt, 's');
        
        if strfind(NOF, 'all')
            FI = 1;
            LI = length(TransfoDir);
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
        
        imNum = TransfoDir(index).name(end-7:end-4);
        
        load([LoadPathPIVRes 'Exp' expName '_PIVRes_' imNum '.mat']);
        load([LoadPathTransfo 'Exp' expName '_transfo_' imNum '.mat']);
        load([LoadPathTransfoFF 'Exp' expName '_transfoFF_' imNum '.mat']);
                
        % Generate transfoH
        SU = imag(transfo.SU(2:end,:));
        transfoH.SU = TransformVelField_decay_hor(transfo.SU(2:end,:), PIVRes, SU);
        
        % Generate TransfoFFH
        SU = imag(transfoFF.SU(2:end,:));
        transfoFFH.SU = TransformVelField_decay_hor(transfoFF.SU(2:end,:), PIVRes, SU);
        transfoFFH.AK = TransformVelField_decay_hor(transfoFF.AK(2:end,:), PIVRes, SU);
        transfoFFH.ORB = TransformVelField_decay_hor(transfoFF.ORB(2:end,:), PIVRes, SU);
        transfoFFH.ORB_acc = TransformVelField_decay_hor(transfoFF.ORB_acc(2:end,:), PIVRes, SU);
        
        FileName = ['Exp' expName '_transfoH_' imNum '.mat'];
        save([SavePathTransfoH FileName], 'transfoH');
        
        FileName = ['Exp' expName '_transfoFFH_' imNum '.mat'];
        save([SavePathTransfoFFH FileName], 'transfoFFH');
        disp(['Pair ' imNum ' done.']);
    
    end
    
end
