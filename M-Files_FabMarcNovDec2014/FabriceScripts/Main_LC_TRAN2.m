TRAN=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Transverse\PIVdt8ms_IRlas1_8hz\']
DIRS=dir(TRAN);
DIRS=DIRS(3:end);

for ii=1:length(DIRS)
exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [TRAN exp_name];
files=dir([load_path '\PIVRaw\PIVCC\*.mat']);
number_of_pair=length(files)/2;

for image_pair_number=0
    [ii  image_pair_number]
%PIV
try
    
load([load_path '\PIVRaw\PIVCC\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = (imrotate(imgPiv,-0.4)); IM_a=IM_a(27:2042,26:2033); %rotation and crop of PIV image
load([load_path '\PIVRaw\PIVCC\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = (imrotate(imgPiv,-0.4)); IM_b=IM_b(27:2042,26:2033);%rotation and crop of PIV image
 
%PIV Surf
load([load_path '\PIVRaw\PIVSURFCC\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurf1 = imgPivsurf;
load([load_path '\PIVRaw\PIVSURFCC\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurf2 = imgPivsurf;

%Surface detection and Creating Masks
U = [224 244; 1904 245; 250 907; 1873 911];
X = [224 244; 1904 244; 224 907; 1904 907];
T = maketform('projective', U, X);
[b1, ~, ~] = imtransform(imgPivsurf1,T,'XYScale',1); %rectify
c1=imresize(b1,[2333 2257]); %Resizing to square
d1=imresize(c1,147.36/94.4); %Resizing to match PIV
dd1=imrotate(d1,-0.55); %rotate
e1=dd1(857:2434,223:3350); %crop to save time
imSurf1 = findSurface_simple_ext_force((medfilt2(e1)), 1); %find surface
Surface_PIV1=imSurf1.surface(518:518+2008);  %in e1 x=518  is the left edge of PIV image
Surface_PIV1=Surface_PIV1-797+433; %still water is Y=797 in e1  and Y=433 in PIV (rotated and cropped)
%imagesc(IM_a); colormap gray, caxis([150 400]); hold; plot(Surface_PIV1,'w')

U = [224 244; 1904 245; 250 907; 1873 911];
X = [224 244; 1904 244; 224 907; 1904 907];
T = maketform('projective', U, X);
[b2, ~, ~] = imtransform(imgPivsurf2,T,'XYScale',1); %rectify
c2=imresize(b2,[2333 2257]); %Resizing to square
d2=imresize(c2,147.36/94.4); %Resizing to match PIV
dd2=imrotate(d2,-0.55); %rotate
e2=dd2(857:2434,223:3350); %crop to save time
imSurf2 = findSurface_simple_ext_force((medfilt2(e2)), 1); %find surface
Surface_PIV2=imSurf2.surface(518:518+2008);  %in e1 x=518  is the left edge of PIV image
Surface_PIV2=Surface_PIV2-797+433; %still water is Y=797 in e1  and Y=433 in PIV (rotated and cropped)


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
compVel.DX=1/14736; %m per pix
compVel.DT=8d-3; % sec per image pair -  DELTA_T= 10 milisec

%post process



%OUTPUT
outfile = [load_path '\PIVMat2\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'compVel', 'imSurf1', 'imSurf2');
disp(['pair ' num2str(image_pair_number) ' done.']);
end
end


end











