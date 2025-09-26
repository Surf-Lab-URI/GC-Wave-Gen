
LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\']
DIRS=dir(LONG);
DIRS=DIRS(3:end);




for ii=1:length(DIRS)

exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = ['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\' exp_name];
files=dir([load_path '\PIVRaw\PIV\*.mat']);
number_of_pair=length(files)/2;

for image_pair_number=0:number_of_pair-1
%PIV
load([load_path '\PIVRaw\PIV\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imgPiv;
load([load_path '\PIVRaw\PIV\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = imgPiv;
 
%PIV Surf
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurf1 = imgPivsurf;
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurf2 = imgPivsurf;

%Surface detection and Creating Masks
d1=imresize(imgPivsurf1,176.77/103.48); %Resizing to match PIV
imSurf1 = findSurface_simple_ext_force((medfilt2(d1)), 1);
Surface_PIV1=imSurf1.surface(724:724+2047); %in this resized surface image, x=724 is the left edge of PIV image
Surface_PIV1=Surface_PIV1-1838+370; %still water is Y=1838 in PIVsurf and Y=370 in PIV

d2=imresize(imgPivsurf2,176.77/103.48);
imSurf2 = findSurface_simple_ext_force((medfilt2(d2)), 1);
Surface_PIV2=imSurf2.surface(724:724+2047);
Surface_PIV2=Surface_PIV2-1838+370;

[h, w] = size(IM_a); %image height and width
mask1=ones(size(IM_a));
mask2=ones(size(IM_a));
for i=1:w
    mask1(1:round(Surface_PIV1(i)),i)=NaN;
    mask2(1:round(Surface_PIV2(i)),i)=NaN;
end

imSurf1.surface=Surface_PIV1;
imSurf1.mask=mask1;
imSurf2.surface=Surface_PIV2;
imSurf2.mask=mask2;

IntrWndw=[256 128 64 32 16 8];
GrdSpc=[128 64 32 16 8 4];

compVel =  computeVelocities_marc_quick_nofilt(IM_a, IM_b, mask1, mask2, IntrWndw, GrdSpc);
compVel.DX=1/17677; %m per pix
compVel.DT=10d-3; % sec per image pair -  DELTA_T= 10 milisec

%post process



%OUTPUT
outfile = [load_path '\PIVMat\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'compVel', 'imSurf1', 'imSurf2');
disp(['pair ' num2str(image_pair_number) ' done.']);
end

end


