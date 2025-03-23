function [imSurf] = CrapperOptimized_FindSurface(img, Step,  S)
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
 
 Sigma1 = 10;  % Sigma used to calculated image derivatives
 Sigma2 = 5;
 Sigma3 = 1; % Sigma used to calculated image derivatives 

Eext1 = ExternalForceImage2D_fab(img,Sigma1).*S;
Eext2 = ExternalForceImage2D_fab(img,Sigma2);
Eext3 = ExternalForceImage2D_fab(img,Sigma3);

Eext1(Eext1<0)=0;
badFrameBool = 0;
warning off
for i=1:size(img,2)

    gv1 = abs(Eext1(:,i));
    
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
    if isempty(locs)
        surface(i) = nan;
        badFrameBool = 1;
    elseif locs + Step > size(img,1) || locs - Step < 1  % bad frame
        badFrameBool = 1;
    else
        gv2 = abs(Eext2(locs-Step:locs+Step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs-Step-1;
    end
    
    clear gv1 gv2 locs

end
temp = surface;
Step = 20;
for i=1:size(img,2)

    locs = temp(i);

    if isempty(locs)
        surface(i) = nan;
        badFrameBool = 1;
    elseif locs + Step > size(img,1) || locs - Step < 1  % bad frame
        badFrameBool = 1;
    else
        gv2 = abs(Eext3(locs-Step:locs+Step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs-Step-1;
    end
    
    clear gv1 gv2 locs

end

% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;

imSurf.surface = surface;
imSurf.badFrameBool = badFrameBool;

end
