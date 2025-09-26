LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Transverse\LIFdt8ms_IRlas1_8hz']
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1:length(DIRS)
exp_name=DIRS(ii).name;
num_of_digits = 4;
load_path = [LONG '\' exp_name];
files=dir([load_path '\IRRaw\*.tif']);
number_of_images=length(files);
fname=files(1).name; fname=fname(1:end-5);

for im=0:number_of_images-1
    
a=imread([load_path '\IRRaw\' fname num2str(im) '.tif']);
a=medfilt2(a);
U = [92 26; 480 26; 47 503; 491 503];
X = [47 26; 491 26; 47 506; 491 503];
T = maketform('projective', U, X);
[b1, ~, ~] = imtransform(a,T,'XYScale',1);
b1=b1(4:513,64:700);
imgIR=imresize(b1,[556 637]);

%NEW correction OLD ONE WAS NOT THE SAME IN X AND Z!!!
imgIR2=imresize(b1,[627 637]);

IR.img=imgIR;
IR.DX=3.880270270270270e-004;
IR.im_number=im+1;

%OUTPUT
outfile = [load_path '\IRMat\' exp_name '_IR_' sprintf(['%0' num2str(num_of_digits) 'd'], im)];
save(outfile, 'IR');
disp(['image ' num2str(im) ' done.']);
end
end


