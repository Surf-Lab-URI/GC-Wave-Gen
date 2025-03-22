%% Bad Pixels PIV Water
IM_tot = 0;
count = 0;
In = 1;
for idx = In:3000
    count = count+1;
    % Indexes for images
    pair_index = (image_index(idx)+1)/2;
    PIV1Dir_temp = PIVAirDir;
    PIV2Dir_temp = PIVWaterDir;
    SurfDir_temp = PIVSurf_LFV_Dir;
    ImageNum_Air1 = PIV1Dir_temp(image_index(idx)).name(max(strfind(PIV1Dir_temp(image_index(idx)).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)).name)-4);
    ImageNum_Air2 = PIV1Dir_temp(image_index(idx)+1).name(max(strfind(PIV1Dir_temp(image_index(idx)+1).name,'_'))+1:length(PIV1Dir_temp(image_index(idx)+1).name)-4);
    ImageNum_Water1 = PIV2Dir_temp(image_index(idx)).name(max(strfind(PIV2Dir_temp(image_index(idx)).name,'_'))+1:length(PIV2Dir_temp(image_index(idx)).name)-4);
    filename = [LoadPath 'PIV Water\' expName '_' runName '_PIV Water_' ImageNum_Water1 '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM_tot = IM_tot+IM1;
    disp(idx)
end
IM_fin = IM_tot/count;
figure;imagesc(IM_fin)

IM__ = IM_fin-medfilt2(IM_fin,[5 5]);
STD = 3; %6
IM_ = IM__;
IM_(IM__>mean(IM__(:))+STD*std(IM__(:))) = NaN;
[r1, c1] = find(isnan(IM_));
r1(c1<10 | c1>4080) = [];
c1(c1<10 | c1>4080) = [];
c1(r1>1090) = [];
r1(r1>1090)=[];
BadPix_Water = [r1,c1];
IM_ = IM_tot/count;
for i = 1:length(BadPix_Water)
    IM_(BadPix_Water(i,1),BadPix_Water(i,2)) = NaN;
end
figure;imagesc(IM_);hold on;plot(BadPix_Water(:,2),BadPix_Water(:,1),'ro')

%% Bad Pixels PIV Air
LoadPath = '\\spray3\d\data\EXPERIMENTS\Dark_after_ExpAW3\Dark_after_ExpAW3_Scene1\RAW\';
FI = 0;
LI = 119;
image_index = FI+1:2:LI;

IM_tot = 0;
count = 0;
In = 1;
for idx = In:LI+1
    count = count+1;
    % Indexes for images
    PairNum = num2str(idx-1,'%.3d');
    filename = [LoadPath 'PIV Air\Dark_after_ExpAW3_Scene1_PIV Air_' PairNum '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM_tot = IM_tot+IM1;
    disp(idx)
end
IM_fin = IM_tot/count;
figure;imagesc(IM_fin)
[r1, c1] = find(IM_fin>20);
BadPix_Air = [r1,c1];

%% Bad Pixels PIVSurf Water
LoadPath = '\\spray3\d\data\EXPERIMENTS\Dark_after_ExpAW3\Dark_after_ExpAW3_Scene1\RAW\';
FI = 0;
LI = 119;
image_index = FI+1:2:LI;

IM_tot = 0;
count = 0;
In = 1;
for idx = In:LI+1
    count = count+1;
    % Indexes for images
    PairNum = num2str(idx-1,'%.3d');
    filename = [LoadPath 'PIVSURF Water\Dark_after_ExpAW3_Scene1_PIVSURF Water_' PairNum '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM_tot = IM_tot+IM1;
    disp(idx)
end
IM_fin = IM_tot/count;
figure;imagesc(IM_fin)
[r1, c1] = find(IM_fin>3);
BadPix_SurfW = [r1,c1];
%% Bad Pixels PIVSurf Air - LFV
LoadPath = '\\spray3\d\data\EXPERIMENTS\Dark_after_ExpAW3\Dark_after_ExpAW3_Scene1\RAW\';
FI = 0;
LI = 59;
image_index = FI+1:2:LI;

IM_tot = 0;
count = 0;
In = 1;
for idx = In:LI+1
    count = count+1;
    % Indexes for images
    PairNum = num2str(idx-1,'%.2d');
    filename = [LoadPath 'PIVSurf Air - LFV\Dark_after_ExpAW3_Scene1_PIVSurf Air - LFV_' PairNum '.raw'];
    [IM1] = (load_Image_IOCoreView_12MP(filename));
    IM_tot = IM_tot+IM1;
    disp(idx)
end
IM_fin = IM_tot/count;
figure;imagesc(IM_fin)
[r1, c1] = find(IM_fin>75);
BadPix_SurfA = [r1,c1];