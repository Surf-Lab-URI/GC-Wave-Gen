% This MATLAB script extracts single-point measurements from WG cameras.

clear
close all
clc

ROOTPath = '\\spray3\d\data\EXPERIMENTS\'; % Spray4

ExpDir = dir([ROOTPath 'Exp*']); % Directory with all the experiments

for i = 1%:length(ExpDir) % Loop on the number of experiments

    ExpAW = ExpDir(i).name(6);
    Acc = ExpDir(i).name(11:14);
    Wind = ExpDir(i).name(17);
    switch Wind
        case '5'
            DeltaT_A = 200d-6;
        case '6'
            DeltaT_A = 120d-6;
        case '7'
            DeltaT_A = 90d-6;
    end


    expName = ['ExpAW' ExpAW '_acc' Acc '_W' Wind 'V' ];

    if strcmp(ExpAW,'6')
        expName = ['ExpAW' ExpAW '_acc' Acc '_W' Wind 'V_LidOpen' ];
    end

    RunDir = dir([ExpDir(i).folder '\' ExpDir(i).name '\Exp*Run*']);
    %%% When we want to analyze the "NoFog" runs, we have to modify RunDir.
    %%% But also, we have to modify the script a little bit to skip PIVAir
    %%% analysis

    for ii = 1%:numel(RunDir) % Loop on the number of Runs

        runName = RunDir(ii).name(20:end);
        if strcmp(ExpAW,'6')
            runName = RunDir(ii).name(28:end);
        end

        % Define Path
        DataPath = [ROOTPath expName '\' expName '_' runName '\' ];
        LoadPath = [DataPath 'RAW\'];
        ResultsPath = [DataPath 'RESULTS_fabio\'];
        SaveWGPath = [ResultsPath 'WG\'];

        if ~exist(SaveWGPath, 'dir')
            mkdir(SaveWGPath);
        end

        %% Parameters
        DX = 1d-2/168.13; % Resolution in m per pix (PIV)
        DT = 1/45; % delta_T between flashes for WG cameras (data rate 45Hz for WG12, 90Hz for WG3)
        WG_up_dir = [LoadPath '\WG12\'];
        WG_dw_dir = [LoadPath '\WG3\'];
        WG12Dir = dir([WG_up_dir '*.raw']);
        WG3Dir = dir([WG_dw_dir '*.raw']);

        %%% NOTE: two constants groups: CST_WG12 (for WG12) and CST_WG34 (for WG34)

        % CST_WG12 %%%%%%%%% CHECK THIS
        if length(expName) == 1
            CST_WG12.ExpName = WG12Dir(1).name(1:6);
        elseif  length(expName) == 2
            CST_WG12.ExpName = WG12Dir(1).name(1:7);
        else
            error('No Experiment name detected')
        end
        CST_WG12.DX = DX;
        CST_WG12.DT = DT;
        CST_WG12.date = WG12Dir(1).date;
        CST_WG12.Total_Time = length(WG12Dir)/43.5;
        CST_WG12.Num_of_WG_images = length(WG12Dir);

        % Define physical parameters
        CST_WG12.AIR_DENSITY = 1.204;     % air density [kg/m3] at 20°C
        CST_WG12.DVISCOSITY = 1.825e-5;   % dynamic viscosity of air [kg/m*s] at 20°C
        CST_WG12.g = 9.81;                 % gravitational acceleration in [m/s2]
        CST_WG12.WATER_DEPTH = 0.7;       % Mean water depth in [m]
        CST_WG12.WATER_DENSITY = 1000;    % in [kg/m3]
        CST_WG12.SURFACE_TENSION = 0.074; % surface tension in [N/m]
        CST_WG12.TOLERANCE = 10e-14;      % numerical tolerance

        % CST_WG34
        CST_WG3.DX = DX;
        CST_WG3.DT = DT/2;
        CST_WG3 = CST_WG12;
        CST_WG3.date = WG3Dir(1).date;
        CST_WG3.Total_Time = length(WG3Dir)*CSTWG3.DT/2;
        CST_WG3.Num_of_WG_images = length(WG3Dir);

        save([ResultsPath 'Movie' expName '_Scene' sceneName '_Parameters.mat'],'CST_WG12','CST_WG34','-append')

        %% Processing frames

        LI = min(length(WG12Dir),length(WG3Dir));
        heightWG2 = nan(1,LI);
        heightWG4 = nan(1,LI);
        for image_index = 2:LI % Main Loop length(FI:LI)

            disp(['WG: processing image ' WG12Dir(image_index).name])

            % Image names
            imagenameWG2 = [WG_up_dir WG12Dir(image_index).name];
            imagenameWG4 = [WG_dw_dir WG3Dir(image_index).name];

            % Extract single-point measurements
            [xWG2(image_index),heightWG2(image_index),xWG4(image_index),heightWG4(image_index)] = extract_WG_signal(imagenameWG2,imagenameWG4);

            disp(heightWG2(image_index))

        end

        WG2.eta=despike_fab(heightWG2);
        WG2.x = xWG2;
        WG4.eta=despike_fab(heightWG4);
        WG4.x = xWG4;

        save([SaveWGPath 'Movie' expName '_Scene' sceneName '_WG.mat'],'WG2','WG4')
    end

end