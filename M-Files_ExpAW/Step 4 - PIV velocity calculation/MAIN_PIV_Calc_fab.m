% This MATLAB script fuses two raw PIV images, PIV 1 and PIV 2 for using as
% an input in PIV calculations.


clear;
clc;

ROOTPath = 'D:\URI_EXP\data\';

 %expName = '4';
 %sceneNameList = '4';
if ~exist('expName', 'var')
    prompt = 'Experiment Name (e.g. 4) = ';
    expName = strsplit(input(prompt, 's'), ','); % The number of exper-
    expName=expName{1};
end

if ~exist('sceneNameList', 'var')
    prompt = 'Scene Name List (e.g. 3) = ';
    sceneNameList = strsplit(input(prompt, 's'), ','); % The number of exper-
    sceneName=sceneNameList{1};
end
         
DataPath=[ROOTPath 'ExpW' expName '\'];
     
     load BF9.dat %list of bad frames for 9s interval
     load BF6.dat %list of bad frames for 6s interval
     load NORM_PIV_SURF.dat %PIV_Surf image normalization;
     load NORM_LFV.dat %LFV image normalization;
        
    LoadPath = [DataPath 'ExpW' expName '_Scene' sceneName '\RAW_DATA'];
    SavePIVPath = [DataPath 'ExpW' expName '_Scene' sceneName '\RESULTS\PIV_Velocities_raw\'];
    SavePIVPath_processed = [DataPath 'ExpW' expName '_Scene' sceneName '\RESULTS\PIV_Velocities_processed\'];
    SaveSurfacesPath = [DataPath 'ExpW' expName '_Scene' sceneName '\RESULTS\Surfaces\'];
    SaveCSTPath = [DataPath 'ExpW' expName '_Scene' sceneName '\RESULTS\'];
    SaveTransfoPath = [DataPath 'ExpW' expName '_Scene' sceneName '\RESULTS\Transformations\'];
    
    if ~exist(SavePIVPath, 'dir')
        mkdir(SavePIVPath);
    end
    if ~exist(SavePIVPath_processed, 'dir')
        mkdir(SavePIVPath_processed);
    end
    if ~exist(SaveSurfacesPath, 'dir')
        mkdir(SaveSurfacesPath);
    end
    if ~exist(SaveCSTPath, 'dir')
        mkdir(SaveCSTPath);
    end
    if ~exist(SaveTransfoPath, 'dir')
        mkdir(SaveTransfoPath);
    end
    
    
   %% Parameters
    IntrWndw = [ 64 32 16 8];
    GrdSpc = [ 32 16 8 4];
    DX=2.54d-2/590.3; %resolution in m per pix.imagesc
    DT=130d-6; %delta_T =130 us between flashes

   PIV1Dir = dir([LoadPath '\F800_upstream\' '*.raw']);
   SURFDir = dir([LoadPath '\PIVSurf\' '*.raw']);
   CST.ExpName=PIV1Dir(1).name(1:12);
   CST.DX=DX;
   CST.DT=DT;
   CST.date=PIV1Dir(1).date;
   CST.Total_Time=length(PIV1Dir)/14.5/2;
   CST.Num_of_PIV_images=length(PIV1Dir);
   CST.Num_of_PIV_pairs=length(PIV1Dir)/2;
   CST.Bad_Frames=BF9;
   CST.IntrWndw=IntrWndw;
   CST.GrdSpc=GrdSpc;

   %% SAVING Experiment Constants & Parameters
   CSTname = ['ExpW' expName '_Scene' sceneName '_Parameters'];
   save([SaveCSTPath CSTname], 'CST')  
   
   %% Fame to process
    % if interupted while processing, look at last filename number
    % that's pair_index; choose FI = pair_index*2-1; LI =
    % length(PIV1Dir)-1;
    % FI = 1; LI = length(PIVDir);
   
        
        disp(['The length of directory is ' num2str(length(PIV1Dir)) '.']);
        prompt = 'PIV Frame Numbers (Must be even e.g. all or 0, 2, 4, etc) = ';
        NOF = input(prompt, 's');
        
        if strfind(NOF, 'all')
            FI = 0;
            LI = length(PIV1Dir)-1;
        else
            NOF = strsplit(NOF, '-');            
            FI = str2double(NOF{1});
            try
                LI = str2double(NOF{2});
             catch
                LI = FI+1;
            end
        end
        

    
    if mod(FI,2)~=0
        FI=FI-1;
    end
       
  %% Processing frames
  
    for image_index = FI+1:2:LI % Main Loop length(FI:LI)
        if ~ismember(image_index-1,BF9)
 %% PIV images
        pair_index=(image_index+1)/2;
        ImageNum_1 = PIV1Dir(image_index).name(max(strfind(PIV1Dir(image_index).name,'_'))+1:length(PIV1Dir(image_index).name)-4);
        ImageNum_2 = PIV1Dir(image_index+1).name(max(strfind(PIV1Dir(image_index+1).name,'_'))+1:length(PIV1Dir(image_index+1).name)-4);
        PairNum = SURFDir(pair_index).name(max(strfind(SURFDir(pair_index).name,'_'))+1:length(SURFDir(pair_index).name)-4);
        %Index is refers to the file number. i.e. first file, second file, etc...
        %Number is the number in the file name. i.e it starts at ZERO
        %for example, the first image is called '0000' index is 1, Num is 0
        
       %first PIV image
            [IM1a] = load_Image_IOCoreView([LoadPath '\F800_upstream\ExpW' expName '_Scene' sceneName' '_F800_upstream_' ImageNum_1 '.raw']);
            [IM1b] = load_Image_IOCoreView([LoadPath '\F801_downstream\ExpW' expName '_Scene' sceneName' '_F801_downstream_' ImageNum_1 '.raw']);
         
       %second PIV image
            [IM2a] = load_Image_IOCoreView([LoadPath '\F800_upstream\ExpW' expName '_Scene' sceneName' '_F800_upstream_' ImageNum_2 '.raw']);
            [IM2b] = load_Image_IOCoreView([LoadPath '\F801_downstream\ExpW' expName '_Scene' sceneName' '_F801_downstream_' ImageNum_2 '.raw']);
            
            [PIV1a_CamAngle_Corrected] = Correct_Angle_PIV_F800_v2(IM1a);
            [PIV1b_CamAngle_Corrected] = Correct_Angle_PIV_F801_v2(IM1b); clear IM1a IM1b
 
            [PIV2a_CamAngle_Corrected] = Correct_Angle_PIV_F800_v2(IM2a);
            [PIV2b_CamAngle_Corrected] = Correct_Angle_PIV_F801_v2(IM2b); clear IM2a IM2b
            %FUSE PIV
            [FusedPIV1] = FusePIV_v2(PIV1a_CamAngle_Corrected,PIV1b_CamAngle_Corrected);
            [FusedPIV2] = FusePIV_v2(PIV2a_CamAngle_Corrected,PIV2b_CamAngle_Corrected);
            clear PIV1a_CamAngle_Corrected PIV1b_CamAngle_Corrected PIV2a_CamAngle_Corrected PIV2b_CamAngle_Corrected
      
            [FusedPIV1] = Pre_process_PIV_Image(FusedPIV1);
            [FusedPIV2] = Pre_process_PIV_Image(FusedPIV2);
            
%% surface image 
try
        [PIVSurf] = load_Image_Jai_PIVSurf([LoadPath '\PIVSurf\ExpW' expName '_Scene' sceneName' '_PIVSurf_' PairNum '.raw']); %load surface image
        [PIVSurf_CR1] = CorrectPIVSurfLensDistortion(PIVSurf); clear PIVSurf
        [PIVSurf_CRR] = CorrectPIVSurf_v2(PIVSurf_CR1);  clear PIVSurf_CR1
        %PIVSurf_CRR contains the PIVsurf image that correspond to FusedPIV - it's
        %matched in size and position.
        
        [BadFramePIVSurf,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface,PIVSurf_PIVMatch] = ExtractSurface_PIVSurf(NORM_PIV_SURF,PIVSurf_CRR,FusedPIV1);
        %Surface detection - 5 outputs,
        %XPIV_PIVSurf_Surface is its x-coordinate in the PIV image (0 is 1st column of PIV)
        %PIVSurf_Surface is the surface detected on the PIVsurf image (a bit longer because the PIVSurf fov is larger than the PIV images)
        %PIVFused_Surface is the surface detected on the PIVsurf image but cropped to match the size and range of the PIVFused images
        %PIVSurf_PIVMatch is the PISsurf image cropped to matches the PIV image; 
        [Mask] = PIVMask(FusedPIV1, PIVFused_Surface);
catch
    continue
end
       
%% Combo surface and phase        
  try
       [LFV] = load_Image_Jai_LFV([LoadPath '\LFV\ExpW' expName '_Scene' sceneName' '_LFV_' PairNum '.raw']);
       LFV_CR1=CorrectLFVLensDistortion(LFV);
       [BadFrameLFVSurf,XPIV_LFV_Surface_PIVMatched,LFV_Surface_PIVMatched,LFV_Surface_Combo_Surface] = ExtractSurface_LFV(NORM_LFV,LFV_CR1,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface);
        %Surface detection - 4 outputs,
        %XPIV_LFV_Surface_PIVMatched is its x-coordinate in the PIV image (0 is 1st column of PIV)
        %LFV_Surface_PIVMatched is the surface detected on the LFV image (way longer)
        %LFV_Surface_Combo_Surface is the surface detected on the LFV image but with the PIVSurf_Surface inserted
       
       FiltLength = 500; MWL=2442; %mean water level
       PixRes.Combo_Surface_X=XPIV_LFV_Surface_PIVMatched;
       PixRes.Combo_Surface_eta = LFV_Surface_Combo_Surface;
       PixRes.Combo_Surface_eta_smth = filtfilt(ones(1,FiltLength)/FiltLength, 1, PixRes.Combo_Surface_eta);
       PixRes.Combo_Surface_eta_smth_phase = angle(hilbert(PixRes.Combo_Surface_eta_smth-MWL));
       
       PixRes.PIVSurf_Surface_eta = PIVSurf_Surface;
       PixRes.PIVSurf_Surface_X  = XPIV_PIVSurf_Surface;
       PixRes.pair_index=pair_index;
       PixRes.ImageNum_1=ImageNum_1;
       PixRes.ImageNum_2=ImageNum_2;
       PixRes.PairNum=PairNum;
       PixRes.ExpName=PIV1Dir(1).name(1:12);
         
       [h, w] = size(FusedPIV1); xpiv = IntrWndw(end)/2:GrdSpc(end):(w - IntrWndw(end)/2);
       i=ismember(PixRes.Combo_Surface_X,xpiv);
       PIVRes.Phase=PixRes.Combo_Surface_eta_smth_phase(i);
       
catch
    continue
end  
       
%% Reject images with bad surface detection
       if (BadFrameLFVSurf==1 || BadFramePIVSurf==1)
           continue
       else
%% PIV Calculations   
        [CompVel] =  ComputeVelocities_Quick_NoFilt_Deform(FusedPIV1, FusedPIV2, Mask, IntrWndw, GrdSpc)
        
        PixRes.PIVFused_Surface=PIVFused_Surface;
        PIVRes.xPIV = CompVel.xPIV; % The x coordinates of center of IntrWndws
        PIVRes.zPIV = CompVel.zPIV; % The y coordinates of center of IntrWndws
        PIVRes.GS = CompVel.GS; % Final grid spacing
        PIVRes.PIVFused_Surface=(PixRes.PIVFused_Surface(PIVRes.xPIV) )/CompVel.GS;
        PIVRes.pair_index=pair_index;
        PIVRes.ImageNum_1=ImageNum_1;
        PIVRes.ImageNum_2=ImageNum_2;
        PIVRes.PairNum=PairNum;
        PIVRes.ExpName=PIV1Dir(1).name(1:12);
        PIVRes.PF_Surface=length(PIVRes.zPIV)-PIVRes.PIVFused_Surface+1; %will need that for transformations; 
        %it's the surface that would be detected on an upside down PIV image

        
        CompVel.Surface=(PixRes.PIVFused_Surface(PIVRes.xPIV) )/CompVel.GS;
        CompVel.pair_index=pair_index;
        CompVel.ImageNum_1=ImageNum_1;
        CompVel.ImageNum_2=ImageNum_2;
        CompVel.PairNum=PairNum;
        CompVel.ExpName=PIV1Dir(1).name(1:12);
        
       
%% SAVING raw PIV data
       PIVFileName = ['ExpW' expName '_Scene' sceneName '_PIV_Velocity_' PairNum];
       save([SavePIVPath PIVFileName], 'CompVel');
%% SAVING raw Surface data      
       PixResName = ['ExpW' expName '_Scene' sceneName '_Surfaces_' PairNum];
       save([SaveSurfacesPath PixResName], 'PixRes', 'PIVRes')   
%% Smoothing PIV data

       [Cartesian] = RemoveOutliers_fab(CompVel);
        Cartesian.PIVSurf=PIVSurf_PIVMatch(CompVel.zPIV,CompVel.xPIV); %save PIVsurf image on PIV grid.
        Cartesian.Surface=(PixRes.PIVFused_Surface(PIVRes.xPIV) )/CompVel.GS;
        Cartesian.pair_index=pair_index;
        Cartesian.ImageNum_1=ImageNum_1;
        Cartesian.ImageNum_2=ImageNum_2;
        Cartesian.PairNum=PairNum;
        Cartesian.ExpName=PIV1Dir(1).name(1:12);
        Cartesian.Phase=PIVRes.Phase;   
        Cartesian.xPIV = CompVel.xPIV; % The x coordinates of center of IntrWndws
        Cartesian.zPIV = CompVel.zPIV; % The y coordinates of center of IntrWndws
        Cartesian.GS = CompVel.GS;
        Cartesian.IW = CompVel.IW;
%% SAVING processed PIV data        
        FileName = ['ExpW' expName '_Scene' sceneName '_CartVel_' PairNum];
        save([SavePIVPath_processed FileName], 'Cartesian', 'PixRes', 'PIVRes', 'CST');

%% SAVING Coordinate transformation data
        [transfo] = GenerateTransfo_Fab(PixRes, PIVRes, CST);
        
        FileName = ['ExpW' expName '_Scene' sceneName '_transfo_' PairNum];
       save([SaveTransfoPath FileName], 'transfo');
            
     
        
%% Screen output
       disp(['Pair ' PairNum ' done.']);
        
       end
       
       end
    
    end
    
   





    
    

