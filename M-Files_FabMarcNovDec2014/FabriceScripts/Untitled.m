tic
TRAN=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Transverse\PIVdt8ms_IRlas1_8hz\ExpLCTB_1_04\IRMat\'];
DIRS=dir(TRAN);
files=DIRS(3:end);

INT=5;


DELX=[];
DELY=[];
DCOR=[];
SUBPIX=[];

s=[10:24:2592-INT];


for i=s
   
filename1=files(i).name;
 [i, i+INT]
filename1=[TRAN filename1];
load(filename1);
IM1=IR.img;

filename2=files(i+INT).name;
filename2=[TRAN filename2];
load(filename2);
IM2=IR.img;



%PIV Cross-corr
bxA = wiener2(IM1(1:end-1,1:end-1));
bxB = wiener2(IM2(1:end-1,1:end-1));
IW=length(bxA);

N=size(bxA);
Nx=N(2);
Ny=N(1);
Hanning2d=0.5*(1-cos(2*pi*(0:Ny-1)'/(Ny-1)))*0.5*(1-cos(2*pi*(0:Nx-1)/(Nx-1)));
                    
bxAmm = bxA-mean(bxA(:));
bxBmm = bxB-mean(bxB(:));
%bxAmm = bxAmm.*Hanning2d;
%bxBmm = bxBmm.*Hanning2d;
fftA = fft2( bxAmm );
fftB = fft2( bxBmm );

fftCorr = fftB .* conj(fftA);
PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
Xcorr  = fftshift(real( ifft2( fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
[PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in Phase correlation
[Xpky, Xpkx]   = find( Xcorr == max(max(Xcorr)) ); %find max in Cross correlation
ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
ldely =  Xpky - IW/2 - 1;            
T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1)); %3 point gaussian interpolation
t = T(:,2);
SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
t = T(2,:);
SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                    
if (isreal(SubpixelY) && isreal(SubpixelX))
delx =  ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
dely =  ldely + SubpixelY;
dcor = max(max(Xcorr));
subpix=1;
else
 delx = ldelx; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
 dely = ldely;
 dcor = max(max(Xcorr));
 subpix=0;
end

DELX=[DELX delx];
DELY=[DELY dely];
DCOR=[DCOR dcor];
SUBPIX=[SUBPIX subpix];
end
DELX=DELX*(IR.DX/(1/43.2*INT));
DELY=DELY*(IR.DX/(1/43.2*INT));
toc





