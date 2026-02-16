LONG=['H:\FabMarcNovDec2014\Data\Longitudinal\LIFdt10ms_IRlas1_8hz'];
DIRS=dir(LONG);
DIRS=DIRS(3:end);
ii=4;
exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = ['H:\FabMarcNovDec2014\Data\Longitudinal\LIFdt10ms_IRlas1_8hz\' exp_name];
files=dir([load_path '\PIVRaw\PIV\*.mat']);
number_of_pair=length(files)/2;


for image_pair_number=120:190
%PIV
load([load_path '\PIVRaw\PIV\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imrotate(imgPiv,-0.13);
IM_a=imgPiv;

imagesc([0 13], [-2.3 10.7], medfilt2(IM_a,[5 5]))
colormap gray; caxis([0 2000])
axis([0 13 -1 5])

 xlabel('cm')
 ylabel('cm')
 
 
 h=gcf;
filestr=['ExpLCTLIF_2_01_LIF_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.tiff'] ;
 print( h, '-dtiff', filestr)
  
end

TRAN=['H:\FabMarcNovDec2014\Data\Transverse\LIFdt8ms_IRlas1_8hz'];
DIRS=dir(TRAN);
DIRS=DIRS(3:end);
ii=2;
exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = ['H:\FabMarcNovDec2014\Data\Transverse\LIFdt8ms_IRlas1_8hz\' exp_name];
files=dir([load_path '\PIVRaw\PIVCC\*.mat']);
number_of_pair=length(files)/2;


for image_pair_number=120:190
%PIV
load([load_path '\PIVRaw\PIVCC\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imrotate(imgPiv,-0.13);
imagesc([0 13], [-2.9 10.1], medfilt2(IM_a,[5 5]))
colormap gray; caxis([0 2000])
axis([0 13 -1 5])

 xlabel('cm')
 ylabel('cm')
 
 
 h=gcf;
filestr=['ExpLCTLIF_2_01_LIF_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.tiff'] ;
 print( h, '-dtiff', filestr)
  
end