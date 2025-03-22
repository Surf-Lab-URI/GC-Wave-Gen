function [imSurf] = FindSurface_Water(img, Step,S)
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
 
 Sigma1 = 25;%20   % Sigma used to calculated image derivatives 
 Sigma2 = Step; % Sigma used to calculated image derivatives 

Eext1 = ExternalForceImage2D_fab(img,Sigma1).*S;
Eext2 = ExternalForceImage2D_fab(img,Sigma2);
Eext1(Eext1<0)=0;
badFrameBool = 0;
warning off
for i=1:size(img,2)

    gv1 = abs(Eext1(:,i));
    
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
    % [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'npeaks',1);
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


%% Post-Correction for sudden jumps
% Despiking jumps
Jump_thr = 20; 
surf2_corr = despike_jumps(surface,Jump_thr,'SURF_WATER');
surface3 = surf2_corr;
% % % 

%% Bad Frame
Length = 101:length(surface)-100; 
SS = surface(Length);
[~,locs1] = findpeaks(abs(diff(SS)), 'minpeakheight', Jump_thr);locs1(locs1<Length(1) | locs1>Length(end)) = [];
SS = surface3(Length);
for ii = 1:length(locs1)
    SS(locs1(ii)-100:locs1(ii)+100) = NaN;
end
idx = isnan(SS);
fit1 = polyfit(Length(~idx),SS(~idx),25);
Surf_smooth = polyval(fit1,Length);
Diff = abs(Surf_smooth-surface3(Length));

Thr_BF = max(1.2*Jump_thr,70);
if max(Diff)>Thr_BF
    badFrameBool = 1;
else
    badFrameBool = 0;
end

%%
% % % surface3 = surface;
imSurf.surface = surface3;
imSurf.badFrameBool = badFrameBool;

end
% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;
