function [ComboPIV] = Combine_PIV(PIV_A,PIV_W, PF_Surface)
% The PIVMask function creates the mask from the PIVFused image and the PIV
% fused surface. Note that the final result for Mask should be flipped left
% to right. Moreover, the PIVFused surface should be on the bottom of image
% not on the top of it.
%
%    PIVFused = pivIn.IM_a;
%    PF_Surface = size(IM_a,1) - pixRes.PF_Surface;

Size1 = size(PIV_W,1); %max(size(PIV_A,1),size(PIV_W,1));
Size2 = size(PIV_W,2); %max(size(PIV_A,2),size(PIV_W,2));

ComboPIV = NaN(Size1,Size2);

for i = 1:size(ComboPIV,2)
    if i <= size(PIV_A,2)
        ComboPIV(1:floor(PF_Surface(i)),i) = PIV_A(1:floor(PF_Surface(i)),i);
    end
    ComboPIV(floor(PF_Surface(i))+1:end,i) = PIV_W(floor(PF_Surface(i))+1:end,i);
end

%For PIVLAB ad lines below
%Mask(Mask==1)=0;
%Mask(isnan(Mask))=1;

end