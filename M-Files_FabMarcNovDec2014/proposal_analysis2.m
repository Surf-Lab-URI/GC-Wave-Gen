clear

TRAN ='D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Transverse\LIFdt8ms_IRlas1_8hz\'

DIRS=dir(TRAN);
DIRS=DIRS(3:end);


for ii=1:length(DIRS)
DD=[];


exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [TRAN exp_name];
files=dir([load_path '\PIVRaw\PIVCC\*.mat']);

number_of_pair=length(files)/2;



for image_pair_number=1:number_of_pair-1
%Load PIV Surf

load([load_path '\PIVMat\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat']); %replace ~ with path


LIFa=compVel.LIFa;
LIFb=compVel.LIFb;
 
interp_LIF=LIFa+LIFb;
interp_LIF=medfilt2(interp_LIF,[8 8]);
% imagesc(medfilt2(interp_LIF(10:end,10:end),[8 8])); title(num2str(image_pair_number))
% pause
DD=[DD nanmean(interp_LIF')']; 
end           
figure
imagesc(DD)
end
