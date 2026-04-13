function surfCropped = CropSurfToPIVDims(surf,verticalCorrection)
%CROPSURFTOPIVDIMS Cropps surface elevation array from size of the surface
%image with small crop to remove lense distortion that has been rescaled to
%PIV coordinates to the size of the PIV image. verticalCorrection offsets
%the surface vertically so that it aligns with the surface in the PIV
%image. If you are using this function on an a surface elevation array that
%has a mean of zero and is not directly being used for PIV, you should set
%this to false, otherwise every value in your array will have a large
%number subtracted from it in addition to the cropping.
surfCropped = surf(655:655+2047);
if verticalCorrection
    surfCropped = surfCropped-1716+287;
end
end

