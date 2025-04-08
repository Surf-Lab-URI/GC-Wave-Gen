function [imSurf] = CrapperOptimized_FindSurface(img, Sigma, Step,  mask)
% This function finds the surface in PIV Surf as well as LFV images. It was
% basically written by Marc Buckley on the August 22 2013, and the original
% file can be found in the following directory:  \MFiles\surface_detection\
% functions.
%
%    This function implements the basic snake segmentation contour. A snake
%    is an active (moving) contour, in which the points are being attracted
%    by edges and other image boundaries. In order to keep  contour smooth,
%    a membrane and thin plate energy are used as contour regularization. 

%
% outputs,
%  Eextern : The energy function described by the image
%

warning off
Eext1 = ExternalForceImage2D_fab(img,Sigma(1)).*mask;
Eext1(Eext1<0)=0;
surface = zeros(1,size(img,2));
badFrameBool = 0;
for i=1:size(img,2)
    gv1 = abs(Eext1(:,i));
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
    if isempty(locs)
        surface(i) = nan;
        badFrameBool = 1;
    elseif locs + Step(1) > size(img,1) || locs - Step(1) < 1  % bad frame
        badFrameBool = 1;
    else
        surface(i) = locs;
    end
    clear gv1 gv2 locs

end

for j = 1:length(Step)
    lastsurf = surface;
    Eext = ExternalForceImage2D_fab(img,Sigma(j+1));
    for i=1:size(img,2)
    
        locs = lastsurf(i);
    
        if isempty(locs)
            surface(i) = nan;
            badFrameBool = 1;
        elseif locs + Step(j) > size(img,1) || locs - Step(j) < 1  % bad frame
            badFrameBool = 1;
        else
            gv2 = abs(Eext(locs-Step(j):locs+Step(j),i));
            [~, s_gv2] = max(gv2);
            surface(i) = s_gv2+locs-Step(j)-1;
        end
        
        clear gv2 locs

    end
end

% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;

imSurf.surface = surface;
imSurf.badFrameBool = badFrameBool;

end
