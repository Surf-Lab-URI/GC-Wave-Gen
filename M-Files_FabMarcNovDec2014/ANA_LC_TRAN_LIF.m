

num_of_digits_PIV = 3;
num_of_digits_IR = 4;

load_path = ['H:\FabMarcNovDec2014\Data\Transverse\LIFdt8ms_IRlas1_8hz\ExpLCTLIF_2_01'];
filesPIV=dir([load_path '\PIVMAT\*.mat']);
filesIR=dir([load_path '\IRMAT\*.mat']);

number_of_PIV=length(filesPIV);
number_of_IR=length(filesIR);

s=0;
for i=1:500
    image_PIV_number=i;
    image_IR_number=i;
    i
%PIV
load([load_path '\IRMat\ExpLCTLIF_2_01_IR_' sprintf(['%0' num2str(num_of_digits_IR) 'd'], image_IR_number) '.mat']); %replace ~ with path
IR.img=double(IR.img); 
 s=s+imgaussian(IR.img,3);
end
s=s/500;
sm=mean(mean(s));
s2=s-sm;


 MASK=ones(511,511);
MASK(1:114,:)=NaN;

for i=120:190%720:2:1150
    image_PIV_number=i;
    image_IR_number=i;
    
%PIV
load([load_path '\PIVMat\ExpLCTLIF_2_01_compVel_' sprintf(['%0' num2str(num_of_digits_PIV) 'd'], image_PIV_number) '.mat']); %replace ~ with path
si=size(compVel.LIFa);
imagesc([0 si(2)*compVel.DX*4],[0 si(2)*compVel.DX*4],compVel.LIFa.*MASK)
  axis([0 si(2)*compVel.DX*4 0.02 0.08])
  caxis([0 3000]); colormap gray
  xlabel('m')
 ylabel('m')


%load([load_path '\IRMat\ExpLCTLIF_2_01_IR_' sprintf(['%0' num2str(num_of_digits_IR) 'd'], image_IR_number) '.mat']); %replace ~ with path
%IR.img=double(IR.img)-s2; si=size(IR.img);
% subplot(2,2,1) 
% imagesc([0 si(2)*IR.DX],[0 si(1)*IR.DX],medfilt2(IR.img))
%  caxis([19.5 20]);
 
 %imagesc([0 si(1)*IR.DX],[0 si(2)*IR.DX], imgaussian(IR.img,3)')
 %caxis([19.7 20]);
 %xlabel('m')
 %ylabel('m')

%  
%  subplot(2,2,2) 
%  s=size(compVel.LIFa);
%  imagesc([0 s(2)*compVel.DX*4],[0 s(2)*compVel.DX*4],compVel.LIFa.*MASK)
%  axis([0 s(2)*compVel.DX*4 0.02 0.08])
%  caxis([0 3000]);
% 
%  subplot(2,2,3) 
%  imagesc([0 s(2)*compVel.DX*4],[0 s(2)*compVel.DX*4],compVel.delx.*MASK*compVel.DX/compVel.DT*100)
%   axis([0 s(2)*compVel.DX*4 0.02 0.08]);
%   caxis([-3 3]);
%  
%   subplot(2,2,4) 
%  imagesc([0 s(2)*compVel.DX*4],[0 s(2)*compVel.DX*4],compVel.dely.*MASK*compVel.DX/compVel.DT*100)
%   axis([0 s(2)*compVel.DX*4 0.02 0.08]);
%   caxis([-3 3]);
%   
 h=gcf;
filestr=['ExpLCTLIF_2_01_LIF_' sprintf(['%0' num2str(num_of_digits_IR) 'd'], image_IR_number) '.tiff'] ;
 print( h, '-dtiff', filestr)
  
end



%post process



%OUTPUT
outfile = [load_path '\PIVMat\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'compVel', 'imSurf1', 'imSurf2');
disp(['pair ' num2str(image_pair_number) ' done.']);
end


end



