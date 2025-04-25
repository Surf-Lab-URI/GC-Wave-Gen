DataPath = '/media/surflab/New Volume/ExpPilot/ExpPilot5/ExpPilot5_Scene2/'; %[ROOTPath 'ExpPilot' expName '/' 'ExpPilot' expName '_Scene' sceneName '/' ];
LoadPath = [DataPath 'RAW/'];
RawDataPath = [DataPath 'RAW/'];
ResultsPath = [DataPath 'RESULTS_Andy/'];

PIVWaterDir = dir([LoadPath 'PIVSURF Water/' '*.raw']); %Same for water


PIV1Dir_temp = PIVWaterDir;

%%
viewing = true;

idx = 0;
while viewing

    idx
    imagename = [PIV1Dir_temp(idx+1).folder '/' PIV1Dir_temp(idx+1).name];
    [IM1] = load_Image_IOCoreView_12MP(imagename);
    figure(50)
    imagesc(IM1,[0,85])
    set(gca,'DataAspectRatio',[1 1 1])
    ip = input("Next Frame?","s")
    ip
    if ip == 'a'
        idx = max(idx-1,0);
    elseif ip == 'd'
        idx = min([size(PIVWaterDir,1)-1,idx+1]);
    else
        try
            ip = uint16(str2double(ip));
            idx = ip;
        end
    end
end


