function [Mask] = PIVWater_Mask(PIV, PF_Surface)
% The PIVMask function creates the mask from the PIVFused image and the PIV
% fused surface. Note that the final result for Mask should be flipped left
% to right. Moreover, the PIVFused surface should be on the bottom of image
% not on the top of it.
% 
%    PIVFused = pivIn.IM_a;
%    PF_Surface = size(IM_a,1) - pixRes.PF_Surface;


Mask = NaN(size(PIV));
    
for i = 1:size(Mask,2)
    Mask(floor(PF_Surface(i))+1:end,i) = 1;
end

%For PIVLAB ad lines below
%Mask(Mask==1)=0;
%Mask(isnan(Mask))=1;

end