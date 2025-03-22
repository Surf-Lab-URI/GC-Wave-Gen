function [imSurf] = FindSurface(img, Step,S,expName)
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
 
 Sigma1 = 20;   % Sigma used to calculated image derivatives 
 Sigma2 = Step; % Sigma used to calculated image derivatives 

Eext1 = ExternalForceImage2D_fab(img,Sigma1).*S;
Eext2 = ExternalForceImage2D_fab(img,Sigma2);
Eext1(Eext1<0)=0;
warning off
Jump_thr = 20; % error if surface jumps by more than <Jump_thr> pixels
HighWindExp = [7:12,19:24,31:36,43:48,49:72];
Thr_BF = max(1.2*Jump_thr,70);
if ismember(str2double(expName),HighWindExp)
    Jump_thr = 20*2;
    Thr_BF = 2*Thr_BF;
end
for i=1:size(img,2)

    gv1 = abs(Eext1(:,i));
    
%     [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'npeaks',1);

    %%%% Modified by Fabio A. for Shoaling waves
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2.5, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
    
    if str2double(expName) > 36
            [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
    end
    
    if str2double(expName) > 48
            [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/3, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
    end
    
%     if str2double(expName) > 24 && str2double(expName)<37
%         if i >4500 && i<6500 %%% if rotated PIVSurf
%             [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
%         end
%     end

    if str2double(expName) <13
        if i >4500 && i<6500 %%% if rotated PIVSurf
            [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
        end
    end

    if str2double(expName) > 12 && str2double(expName)<37
%         if i >4500 && i<6500 %%% if rotated PIVSurf
            [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'MinPeakProminence', 1e-2); %, 'MinPeakDistance', 100, 'npeaks',1);
%         end
    end
    
    %%%%
    if isempty(locs)
        surface(i) = nan;
        surface2(i) = nan;
        badFrameBool = 1;
    elseif locs(1) + Step > size(img,1) || locs(1) - Step < 1  % bad frame
        badFrameBool = 1;
    else
        gv2 = abs(Eext2(locs(1)-Step:locs(1)+Step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs(1)-Step-1;
        surface2(i) = locs(1);
    end
    
    % If there's a jump, recalculate local peaks and find peak that minimizes 
    % distance with previous point on surface. If distance is still more
    % than expected, interpolate surface based on previous trend
    
    clear gv1 gv2 locs

end

%% Post-Correction for sudden jumps
% Despiking jumps in Shoaling Waves
surf2_corr = despike_jumps(surface,Jump_thr,'PIVSURF');
surface3 = surf2_corr;

%% Last attempt to recover frame before marking as bad frame
Length = 1501:9300;
SS = surface(Length);
[~,locs1] = findpeaks(abs(diff(SS)), 'minpeakheight', Jump_thr);locs1(locs1<Length(1) | locs1>Length(end)) = [];
SS = surface3(Length);
for ii = 1:length(locs1)
    SS(locs1(ii)-100:locs1(ii)+100) = NaN;
end
idx = isnan(SS);
fit1 = polyfit(Length(~idx),SS(~idx),50);
if str2double(expName)>36 && str2double(expName)<49
    fit1 = polyfit(Length(~idx),SS(~idx),20);
end
Surf_smooth = polyval(fit1,Length);
Surf_smooth = [surface3(1:Length(1)+500-1) Surf_smooth(501:end-500) surface3(Length(end)-500+1:end)];
Diff = abs(Surf_smooth-surface3);
Int = 1:length(Diff);
Int(Diff<Jump_thr) = [];
Int(isnan(surface3(Int))) = [];
IntP = Int+100;
IntM = Int-100;
Int = unique([Int,IntP,IntM]);
%
for i = 1:length(Int)
    gv1 = abs(Eext1(:,Int(i)));
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/6, 'MinPeakProminence', 1e-2); %, 'npeaks',1);
    if ~isempty(locs)
        [~,I] = min(abs(locs-Surf_smooth(Int(i))));
        gv2 = abs(Eext2(locs(I)-Step:locs(I)+Step,Int(i)));
        [~, s_gv2] = max(gv2);
        surface3(Int(i)) = s_gv2+locs(I)-Step-1;
    end
end
surf2_corr = despike_jumps(surface3,Jump_thr,'PIVSURF');
surface3 = surf2_corr;

%% Bad Frame
Length = 1501:9300;
SS = surface(Length);
[~,locs1] = findpeaks(abs(diff(SS)), 'minpeakheight', Jump_thr);locs1(locs1<Length(1) | locs1>Length(end)) = [];
SS = surface3(Length);
for ii = 1:length(locs1)
    SS(locs1(ii)-100:locs1(ii)+100) = NaN;
end
idx = isnan(SS);
fit1 = polyfit(Length(~idx),SS(~idx),25);
if str2double(expName) > 36
    fit1 = polyfit(Length(~idx),SS(~idx),10);
end
Surf_smooth = polyval(fit1,Length);
Diff = abs(Surf_smooth-surface3(Length));

if max(Diff)>Thr_BF
    badFrameBool = 1;
else
    badFrameBool = 0;
end

%%
imSurf.surface = surface3;
imSurf.badFrameBool = badFrameBool;

end
