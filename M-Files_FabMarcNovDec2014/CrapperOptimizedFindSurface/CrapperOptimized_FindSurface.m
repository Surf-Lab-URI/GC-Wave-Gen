function [imSurf] = CrapperOptimized_FindSurface(img, Sigma, Step,  mask, slopeDiffThreshold)
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
% Eext1 = ExternalForceImage2D_fab_ExpAW(img,Sigma(1)).*mask;
Eext1 = -ExternalForceImage2D(img,0,1,0,Sigma(1)).*mask;
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
    %Eext = ExternalForceImage2D_fab_ExpAW(img,Sigma(j+1));
    Eext = ExternalForceImage2D(img,0,0.8,0.2,Sigma(j+1));
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

    % Check for weird jumps introduced by the latest step. If the slope
    % of the surface at certain location is way higher than at the 
    % previous sigma, don't use this or subsequent sigma refinements.

    if any(abs(diff(surface)-diff(lastsurf)) > slopeDiffThreshold)
        surface = lastsurf;
        disp(['Exceeded slopeDiffThreshold at Sigma ',int2str(Sigma(j+1))])
        break
    end
end

% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;

% inn = ~isnan(surface);
% i1 = (1:numel(surface)).';
% pp = interp1(i1(inn),surface(inn),'linear','pp');
% surface_s = fnval(pp,linspace(i1(1),i1(end),length(surface)))';

imSurf.surface_raw = surface;
% imSurf.surface = surface_s;
imSurf.badFrameBool = badFrameBool;
end
