

num_of_digits_PIV = 3;
load_path = ['H:\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\ExpLCL_2_02'];
filesPIV=dir([load_path '\PIVMat\*.mat']);
number_of_PIV=length(filesPIV);

U=[];
OM=[];
SH=[];
U2=[];

for i=120:190
    image_PIV_number=i;
      
%PIV
load([load_path '\PIVMat\ExpLCL_2_02_compVel_' sprintf(['%0' num2str(num_of_digits_PIV) 'd'], image_PIV_number) '.mat']); %replace ~ with path

u=compVel.delta_x;
w=compVel.delta_z;
c=compVel.dcor;
mask=compVel.mask;
u(c<0.4)=NaN; w(c<0.4)=NaN;

U2=[U2   nanmean((u.*DX/DT).^2+(w.*DX/DT).^2,2)];
u=smoothn(u,0.4,'robust');
w=smoothn(w,0.4,'robust');
x=compVel.xPIV;
z=compVel.zPIV;
[dudx, dudz] = csapsDiff(u, 0.001, x,z);
[dwdx, dwdz] = csapsDiff(w, 0.001, x,z);

GS=compVel.GS;
DX=compVel.DX;
DT=compVel.DT;

om2=sqrt(((dudz-dwdx).*mask).^2)/DT;
%imagesc([0 s(2)*DX*GS*100],[-2 s(1)*DX*GS*100-2],om2)
%caxis([0 40]);
s=size(u);
sh2=sqrt((((dudx-dwdz).*mask).^2+((dudz+dwdx).*mask).^2))/DT;
imagesc([0 s(2)*DX*GS*100],[-2 s(1)*DX*GS*100-2],sh2)
caxis([0 50]);

u.*mask.*DX/DT;
%s=size(u);
%imagesc([0 s(2)*DX*GS*100],[-2 s(1)*DX*GS*100-2],u.*mask.*DX/DT*100)
%caxis([-1 20]);
 
 U=[U   nanmean(u,2)];
 OM=[OM   nanmean(om2,2)];
 SH=[OM   nanmean(sh2,2)];
  
     xlabel('cm')
 ylabel('cm')
 
 
 h=gcf;
filestr=['ExpLCTLIF_2_01_LIF_' sprintf(['%0' num2str(num_of_digits_PIV) 'd'], i) '.tiff'] ;
 print( h, '-dtiff', filestr)
 
end






