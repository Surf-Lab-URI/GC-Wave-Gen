% Note from Andy: This script (not Main_IR_LON) was the one used for
% generating the .mat files in the IRMat directory and is believed to be
% correct.

%LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\LIFdt10ms_IRlas1_8hz\']
LONG = ['/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/'];
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1%:length(DIRS)
exp_name=DIRS(ii).name;
num_of_digits = 4;
load_path = [LONG  exp_name];
 files=dir([load_path '/IRRaw/*.tif']);
number_of_images=length(files);
fname=files(1).name; fname=fname(1:end-5);

for im=0%:number_of_images-1
a=imread([load_path '/IRRaw/' fname num2str(im) '.tif']);
a=medfilt2(double(a));
U = [94 25; 480 22; 51 500; 495 501];
X = [51 22; 495 22; 51 501; 495 501];
T = maketform('projective', U, X); 
[b1, ~, ~] = imtransform(a,T,'XYScale',1); % Rectification for IR
b1=b1(7:518,64:698);  %croping
imgIR=imresize(b1,[635 635]); % making calibration same in x and z

IR.img=imgIR;
IR.DX=3.87e-004;
IR.im_number=im+1;

%OUTPUT
outfile = [load_path '\IRMat\' exp_name '_IR_' sprintf(['%0' num2str(num_of_digits) 'd'], im)];
save(outfile, 'IR');
disp(['image ' num2str(im) ' done.']);
end
end


LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\LIFdt10ms_IRlas1_8hz\']
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1:length(DIRS)
exp_name=DIRS(ii).name;
num_of_digits = 4;
load_path = [LONG  exp_name];
files=dir([load_path '\IRRaw\*.tif']);
number_of_images=length(files);
fname=files(1).name; fname=fname(1:end-5);

for im=0:number_of_images-1
    
a=imread([load_path '\IRRaw\' fname num2str(im) '.tif']);
a=medfilt2(double(a));
U = [94 25; 480 22; 51 500; 495 501];
X = [51 22; 495 22; 51 501; 495 501];
T = maketform('projective', U, X); 
[b1, ~, ~] = imtransform(a,T,'XYScale',1); % Rectification for IR
b1=b1(7:518,64:698);  %croping
imgIR=imresize(b1,[635 635]); % making calibration same in x and z

IR.img=imgIR;
IR.DX=3.87e-004;
IR.im_number=im+1;

%OUTPUT
outfile = [load_path '\IRMat\' exp_name '_IR_' sprintf(['%0' num2str(num_of_digits) 'd'], im)];
save(outfile, 'IR');
disp(['image ' num2str(im) ' done.']);
end
end

