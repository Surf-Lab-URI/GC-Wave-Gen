function surfCropped = CropSurfToPIVDims(surf,verticalCorrection)
%CROPSURFTOPIVDIMS Cropps surface elevation array from size of the surface
%image with small crop to remove lense distortion that has been rescaled to
%PIV coordinates to the size of the PIV image
surfCropped = surf(655:655+2047);
if verticalCorrection
    surfCropped = surfCropped-1716+287;
end
end

