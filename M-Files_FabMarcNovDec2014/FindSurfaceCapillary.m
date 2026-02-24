function Surf = FindSurfaceCapillary(path,NameValueArgs)
    arguments
        path
        NameValueArgs.findMask = false;
        NameValueArgs.returnImgs = true;
    end
%FINDSURFACECAPILLARY Finds the surface in a surface image
%   Takes as input a path and to the surface image. If you want want it to
%   spit out a mask with dimensions of the PIV image that masks off the
%   air, you can specify findMask = true as an arguement. The output struct
%   has the following images:
%       ImgScaledCroppedToPIV:  Surface image scaled and cropped to have 
%                               the same dimensions as the PIV image
%       ImgScaledToPIVSmallCrop: Surface image scaled to have the same
%                           resolution as the PIV image, but not the same
%                           dimensions. This image captures more of the
%                           surface than the PIV image, although it is
%                           cropped slightly to eliminate artifacts of lens
%                           distortion correction. Surface detection
%                           is run on this image with output 
%                           surfaceSurfImgScaled. The surface is then
%                           cropped to produce surfacePIVImg, which has the
%                           same dimensions as the PIV image. 
%        ImgScaledToPIV:    Surface images scaled to have the same
%                           resolution as the PIV image but not cropped at 
%                           all.
%
%
%   There are also the following surfaces:
%       surfaceSurfImgScaled:   the filtered surface output by the
%                               surface detection function run on 
%                               ImgScaledToPIVSmallCrop.
%       surfacePIVImg:          surfaceSurfImgScaled cropped down to the
%                               size of the PIV image
%       surface_raw:            Unfiltered surface from surface detection,
%                               has same dimensions as
%                               surfaceSurfImgScaled.

load(path,"imgPivsurf"); %replace ~ with path

[scaledImg,scaledImgSmallCrop,scaledCroppedImg] = SurfImgToPIVDims(imgPivsurf);

surfSigmas = [50 40 30 20 15];
surfSteps = [50 40 30 5];
% SurfMask = 1;
slopeDiffThreshold = 5;
Surf = CrapperOptimized_FindSurface(scaledImgSmallCrop, surfSigmas, surfSteps, 1, slopeDiffThreshold);
Surf.surfaceSurfImgScaled = FiltSurf(Surf.surface_raw,200);

Surf.surfacePIVImg = CropSurfToPIVDims(Surf.surfaceSurfImgScaled);

if NameValueArgs.findMask
    w = 2048;
    mask=ones(w);

    for i=1:w
        mask(1:round(Surf.surfacePIVImg(i)),i)=NaN;
    end

    warning off
    Surf.mask=mask;
end

if NameValueArgs.returnImgs
    Surf.ImgScaledCroppedToPIV = scaledCroppedImg;
    Surf.ImgScaledToPIVSmallCrop = scaledImgSmallCrop;
    Surf.ImgScaledToPIV = scaledImg;
end

end

