clear

LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\LIFdt10ms_IRlas1_8hz']
DIRS=dir(LONG);
DIRS=DIRS(3:end);


for ii=1:length(DIRS);  
exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [LONG,'\' exp_name];
files=dir([load_path '\PIVRaw\PIVSURF\*.mat']);
number_of_pair=length(files)/2;
DX=1/17677; %m per pix for PIV image resolution   

A1=[];
A2=[];
for image_pair_number=0:number_of_pair-1
 [ii image_pair_number]
%PIV LFV to get the surface
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurf1 = imgPivsurf;
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurf2 = imgPivsurf;

%Detecting the surface to calculate wavenumber
d1=imresize(imgPivsurf1,176.77/103.48); %Resizing to match PIV
imSurf1 = findSurface_simple_ext_force((medfilt2(d1)), 1); 
surface_profile1=imSurf1.surface*DX; surface_x1=[1:length(imSurf1.surface)]*DX; %scaling  surface profile in m
d2=imresize(imgPivsurf2,176.77/103.48);
imSurf2 = findSurface_simple_ext_force((medfilt2(d2)), 1);
surface_profile2=imSurf2.surface*DX; surface_x2=[1:length(imSurf2.surface)]*DX; %scaling  surface profile in m

outfile = [load_path '\PIVRaw\EXTRACTED_SURFACES\' exp_name '_Surface_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'surface_profile1', 'surface_x1', 'surface_profile2', 'surface_x2');
disp(['pair ' num2str(image_pair_number) ' done.']);

A1=[A1 surface_profile1];
A2=[A2 surface_profile2];

end

save([load_path '\PIVRaw\EXTRACTED_SURFACES\A'], 'A1', 'A2')

end



        
        




