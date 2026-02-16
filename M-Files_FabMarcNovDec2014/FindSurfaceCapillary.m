function imSurf = FindSurfaceCapillary(path,NameValueArgs)
    arguments
        path
        NameValueArgs.maskDims = NaN; 
    end
%FINDSURFACECAPILLARY Summary of this function goes here
%   Detailed explanation goes here

load(path,"imgPivsurf"); %replace ~ with path

%Surface detection and Creating Masks
U1 = [147 49;2024 57; 1995 1004; 161 999];
X1 = [147 49; 2024 49; 2024 1004; 147 1004];
T1 = fitgeotrans(U1,X1,'projective');
d=imwarp(imgPivsurf,T1,'cubic');

d=imresize(d,176.9769/105.5880);%Resizing to match PIV

s=d(30:3525,755:755+2047); %cropping

surfSigmas = [50 40 30 20 15];
surfSteps = [50 40 30 5];
% SurfMask = 1;
slopeDiffThreshold = 5;
imSurf = CrapperOptimized_FindSurface(s, surfSigmas, surfSteps, 1, slopeDiffThreshold);
imSurf.surface = FiltSurf(imSurf.surface_raw,200);

imSurf.surfacePreOffset = imSurf.surface;
imSurf.surface=imSurf.surface-1716+287;

if ~isnan(NameValueArgs.maskDims)
    w = NameValueArgs.maskDims(2);
    mask=ones(NameValueArgs.maskDims);

    for i=1:w
        mask(1:round(imSurf.surface(i)),i)=NaN;
    end

    warning off
    imSurf.mask=mask;
end

imSurf.s = s;
imSurf.d = d;

end

